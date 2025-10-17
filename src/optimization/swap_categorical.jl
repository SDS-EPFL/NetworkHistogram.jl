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
                counts[g1, g2] += 1
                realized[g1, g2][e] += 1
            else
                break
            end
        end
    end

    for g2 in 1:n_groups, g1 in g2:n_groups
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

    # Pre-allocate workspace with copies of current state
    w = WorkspaceDiscreteSwap(
        SymArray(n_groups, 0.0),
        SymArray(n_groups, 0),
        SymArray(n_groups, zeros(Float64, M)),
        SymArray(n_groups, zeros(Float64, M))
    )

    # Create assignment first
    assignment = Assignment(
        node_labels, edge_list, dists, θ, log_likelihood_per_group, w)

    # Now copy the actual workspace data into w
    for g2 in 1:n_groups, g1 in g2:n_groups
        w.log_likelihood_per_group[g1, g2] = log_likelihood_per_group[g1, g2]
        w.counts[g1, g2] = counts[g1, g2]
        copyto!(w.realized[g1, g2], realized[g1, g2])
        copyto!(w.estimated[g1, g2], estimated[g1, g2])
    end

    return assignment
end

function make_workspace(a::Assignment{E, Dist{D},
        F, W}) where {E, F, D <: Cat, W}
    # Pre-allocate workspace instead of deepcopy
    k = number_groups(a)
    m = num_categories(unwrap(a.θ[1, 1]))

    log_ll = SymArray(k, 0.0)
    counts = SymArray(k, 0)
    realized = SymArray(k, zeros(Float64, m))
    estimated = SymArray(k, zeros(Float64, m))

    return WorkspaceDiscreteSwap(log_ll, counts, realized, estimated)
end

function copy_categorical_workspace!(
        dest::WorkspaceDiscreteSwap, src_assignment::Assignment)
    # In-place copy without allocation
    copy!(dest.log_likelihood_per_group, src_assignment.log_likelihood)

    # Copy counts (scalars)
    copy!(dest.counts, src_assignment.additional_workspace.counts)

    # Copy vector-valued SymArrays element by element
    # Use sparse matrix iteration instead of .d dictionary
    k = size(dest.realized, 1)
    copy_with_array!(dest.realized, src_assignment.additional_workspace.realized)
    copy_with_array!(dest.estimated, src_assignment.additional_workspace.estimated)

    # @inbounds for j in 1:k, i in 1:j
    #     copyto!(dest.realized[i, j], src_ws.realized[i, j])
    # end

    # @inbounds for j in 1:k, i in 1:j
    #     copyto!(dest.estimated[i, j], src_ws.estimated[i, j])
    # end
end

function make_swap_workspace!(ws::WorkspaceDiscreteSwap, a::Assignment)
    # Use in-place copy instead of deepcopy
    copy_categorical_workspace!(ws, a)
end

function revert_swap_workspace!(a::Assignment, ws::WorkspaceDiscreteSwap)
    # Use in-place copy instead of deepcopy
    copy!(a.log_likelihood, ws.log_likelihood_per_group)

    copy!(a.additional_workspace.log_likelihood_per_group, ws.log_likelihood_per_group)
    copy!(a.additional_workspace.counts, ws.counts)

    # Copy vector-valued SymArrays element by element
    # Use sparse matrix iteration instead of .d dictionary
    k = size(ws.realized, 1)
    copy_with_array!(a.additional_workspace.realized, ws.realized)
    copy_with_array!(a.additional_workspace.estimated, ws.estimated)
    # @inbounds for j in 1:k, i in 1:j
    #     copyto!(as.realized[i, j], ws.realized[i, j])
    # end

    # @inbounds for j in 1:k, i in 1:j
    #     copyto!(as.estimated[i, j], ws.estimated[i, j])
    # end
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
        as.additional_workspace.counts[gu, g_inter] -= 1
        as.additional_workspace.realized[gu, g_inter][e] -= 1
        as.additional_workspace.counts[gv, g_inter] += 1
        as.additional_workspace.realized[gv, g_inter][e] += 1
    end
    for (node, e) in iterate_neighbors(as.edges, v)
        if node == u
            continue
        end
        g_inter = as.node_labels[node]
        as.additional_workspace.counts[gv, g_inter] -= 1
        as.additional_workspace.realized[gv, g_inter][e] -= 1
        as.additional_workspace.counts[gu, g_inter] += 1
        as.additional_workspace.realized[gu, g_inter][e] += 1
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
    copy!(as.log_likelihood, as.additional_workspace.log_likelihood_per_group)
    # as.log_likelihood = deepcopy(as.additional_workspace.log_likelihood_per_group)
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
