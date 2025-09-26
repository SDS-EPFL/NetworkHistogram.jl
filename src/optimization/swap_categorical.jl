mutable struct WorkspaceDiscreteSwap{
    D, C <: SymArray, R <: SymArray,
    R2 <: SymArray, S <: SymArray{D},
    L <: SymArray}
    θ::S
    log_likelihood_per_group::L
    counts::C
    realized::R
    estimated::R2
end

struct Cat{M, V <: AbstractVector{<:Real}}
    p::V
    function Cat(p::AbstractVector{<:Real})
        new{Val{length(p)}, typeof(p)}(p / sum(p))
    end
end

function Base.show(io::IO, c::Cat)
    print(io, "Cat($(c.p))")
end

num_categories(::Type{Cat{Val{M}, V}}) where {M, V} = M
num_categories(::Cat{Val{M}, V}) where {M, V} = M
zero(c::Cat{Val{M}, V}) where {M, V} = Cat(ones(eltype(V), M))
distance(c1::Cat{M, V}, c2::Cat{M, V}) where {M, V} = sum(abs.(c1.p .- c2.p))
eltype(::Cat{M, V}) where {M, V} = Int
params(c::Cat{M, V}) where {M, V} = (c.p,)
logpdf(c::Cat, x::Int) = log(c.p[x])
function fit(::Cat{Val{M}, V}, x::AbstractVector{Int}) where {M, V}
    p_est = zeros(eltype(V), M)
    for xi in x
        p_est[xi] += 1
    end
    return Cat{Val{M}, V}(p_est ./ length(x))
end

function sample(c::Cat{Val{M}, V}) where {M, V}
    return findfirst(x -> x >= rand(), cumsum(c.p))
end

function fit(::Cat{Val{M}, V}, x::Int) where {M, V}
    p_est = zeros(eltype(V), M)
    p_est[x] = 1.0
    return Cat(p_est)
end

function set_params!(c::Dist{Cat{M, V}}, p::V) where {M, V}
    set_params!(c.dist, p)
end
function set_params!(c::Cat{M, V}, p::V) where {M, V}
    c.p .= p
end

function Assignment(
        node_labels, edge_list::EdgeList{E},
        dist::Dist{D}) where {E, D <: Cat}
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
            set_params!(θ[g1, g2], estimated[g1, g2])
            for m in 1:num_categories(unwrap(dist))
                if realized[g1, g2][m] > 0
                    log_likelihood_per_group[g1, g2] += realized[g1, g2][m] *
                                                        logpdf(θ[g1, g2], m)
                end
            end
        end
    end
    w = WorkspaceDiscreteSwap(deepcopy(θ), deepcopy(log_likelihood_per_group),
        counts, deepcopy(realized), deepcopy(estimated))
    return Assignment(
        node_labels, edge_list, dists, θ, log_likelihood_per_group, w)
end

function make_workspace(a::Assignment{E, Dist{D},
        F, W}) where {E, F, D <: Cat, W}
    return deepcopy(a.additional_workspace)
end

function make_swap_workspace!(ws::WorkspaceDiscreteSwap, a::Assignment)
    ws.θ = deepcopy(a.θ)
    ws.log_likelihood_per_group = deepcopy(a.log_likelihood)
    ws.realized = deepcopy(a.additional_workspace.realized)
    ws.estimated = deepcopy(a.additional_workspace.estimated)
end

function revert_swap_workspace!(a::Assignment, ws::WorkspaceDiscreteSwap)
    a.θ = deepcopy(ws.θ)
    a.log_likelihood = deepcopy(ws.log_likelihood_per_group)
    as = a.additional_workspace
    as.θ = deepcopy(ws.θ)
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

    for g2 in 1:n_groups
        for g1 in g2:n_groups
            set_params!(as.additional_workspace.θ[g1, g2],
                as.additional_workspace.estimated[g1, g2])
            as.additional_workspace.log_likelihood_per_group[g1, g2] = _fast_ll(
                as.additional_workspace.estimated[g1, g2], as.additional_workspace.realized[
                    g1, g2],
                as.additional_workspace.counts[g1, g2])
        end
    end

    as.θ = deepcopy(as.additional_workspace.θ)
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

function _fast_ll(
        p::AbstractVector, r::AbstractVector, c::Real)
    ll = zero(eltype(p))
    if c > 0
        @inbounds for m in eachindex(p)
            if r[m] > 0
                ll += r[m] * log(p[m])
            end
        end
    end
    return ll
end
