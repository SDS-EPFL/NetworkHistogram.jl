mutable struct WorkspaceDiscreteSwap{C <: SymArray, R <: SymArray,
    R2 <: SymArray, L <: SymArray}
    log_likelihood_per_group::L
    counts::C
    realized::R
    estimated::R2
end

function Assignment(
        node_labels, edge_list::EdgeList{E},
        dist::Dist{Cat{M, T}}) where {E, M, T}
    n_groups = length(unique(node_labels))
    n_nodes = length(node_labels)
    dists = fit(dist, edge_list)
    realized = SymArray(n_groups, zeros(Float64, num_categories(unwrap(dist))))
    estimated = SymArray(n_groups, zeros(Float64, num_categories(unwrap(dist))))
    counts = SymArray(n_groups, 0)

    for u in 1:n_nodes
        g1 = node_labels[u]
        for (v, e) in iterate_neighbors(edge_list, u)
            g2 = node_labels[v]
            if v < u
                counts[minmax(g1, g2)...] += 1
                realized[minmax(g1, g2)...][e] += 1
            else
                break
            end
        end
    end

    for g2 in 1:n_groups, g1 in g2:n_groups
        counts[g1, g2] = counts[minmax(g1, g2)...]
        realized[g1, g2] = realized[minmax(g1, g2)...]
        _fast_normalization!(
            estimated[g1, g2], realized[g1, g2], counts[g1, g2])
    end

    θ = SymArray(n_groups, zero(dist))
    log_likelihood_per_group = SymArray(n_groups, 0.0)
    for g2 in 1:n_groups
        for g1 in g2:n_groups
            θ[g1, g2] = Dist(Cat(SVector{M}(estimated[g1, g2])))
            log_likelihood_per_group[g1, g2] = logpdf_cat(
                estimated[g1, g2], realized[g1, g2])
        end
    end
    w = WorkspaceDiscreteSwap(deepcopy(log_likelihood_per_group),
        counts, deepcopy(realized), deepcopy(estimated))
    return Assignment(
        node_labels, edge_list, dists, θ, log_likelihood_per_group, w)
end

function make_workspace(a::Assignment{E, Dist{D},
        F, W}) where {E, F, D <: Cat, W}
    return deepcopy(a.additional_workspace)
end

function make_swap_workspace!(ws::WorkspaceDiscreteSwap, a::Assignment)
    ws.log_likelihood_per_group = deepcopy(a.log_likelihood)
    ws.realized = deepcopy(a.additional_workspace.realized)
    ws.estimated = deepcopy(a.additional_workspace.estimated)
end

function revert_swap_workspace!(a::Assignment, ws::WorkspaceDiscreteSwap)
    a.log_likelihood = deepcopy(ws.log_likelihood_per_group)
    as = a.additional_workspace
    as.log_likelihood_per_group = deepcopy(ws.log_likelihood_per_group)
    as.realized = deepcopy(ws.realized)
    as.estimated = deepcopy(ws.estimated)
end

function apply_swap!(as::Assignment, s::Swap{<:WorkspaceDiscreteSwap})
    u, v = s.u, s.v
    n_groups = number_groups(as)
    gu = as.node_labels[u]
    gv = as.node_labels[v]
    for (node, e) in iterate_neighbors(as.edges, u)
        if node == v
            continue
        end
        g_inter = as.node_labels[node]
        as.additional_workspace.counts[minmax(gu, g_inter)...] -= 1
        as.additional_workspace.realized[minmax(gu, g_inter)...][e] -= 1
        as.additional_workspace.counts[minmax(gv, g_inter)...] += 1
        as.additional_workspace.realized[minmax(gv, g_inter)...][e] += 1
    end
    for (node, e) in iterate_neighbors(as.edges, v)
        if node == u
            continue
        end
        g_inter = as.node_labels[node]
        as.additional_workspace.counts[minmax(gv, g_inter)...] -= 1
        as.additional_workspace.realized[minmax(gv, g_inter)...][e] -= 1
        as.additional_workspace.counts[minmax(gu, g_inter)...] += 1
        as.additional_workspace.realized[minmax(gu, g_inter)...][e] += 1
    end
    _fast_normalization!.(as.additional_workspace.estimated,
        as.additional_workspace.realized, as.additional_workspace.counts)
    swap_node_labels!(as, u, v)
    m = size(as.additional_workspace.estimated[1, 1], 1)
    for g2 in 1:n_groups
        for g1 in g2:n_groups
            as.θ[g1, g2] = Dist(Cat(SVector{m}(as.additional_workspace.estimated[g1, g2])))
            # set_params!(as.additional_workspace.θ[g1, g2],
            #     as.additional_workspace.estimated[g1, g2])
            as.additional_workspace.log_likelihood_per_group[g1, g2] = logpdf_cat(
                as.additional_workspace.estimated[g1, g2], as.additional_workspace.realized[
                    g1, g2])
        end
    end

    as.log_likelihood = deepcopy(as.additional_workspace.log_likelihood_per_group)
end

function _fast_normalization!(p::AbstractVector, r::AbstractVector, c::Real)
    if c > 0
        @inbounds for m in eachindex(p)
            p[m] = r[m] / c
        end
    else
        fill!(p, 0.0)
    end
end
