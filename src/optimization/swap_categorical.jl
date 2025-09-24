mutable struct WorkspaceDiscreteSwap{
    D, C <: AbstractMatrix, R <: AbstractMatrix,
    R2 <: AbstractMatrix, S <: AbstractMatrix{D},
    L <: AbstractMatrix}
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

num_categories(::Type{Cat{Val{M}, V}}) where {M, V} = M
num_categories(::Cat{Val{M}, V}) where {M, V} = M
zero(c::Cat{Val{M}, V}) where {M, V} = Cat(ones(eltype(V), M))
distance(c1::Cat{M, V}, c2::Cat{M, V}) where {M, V} = sum(abs.(c1.p .- c2.p))
eltype(::Cat{M, V}) where {M, V} = Int
params(c::Cat) = (c.p,)
logpdf(c::Cat, x::Int) = log(c.p[x])

function set_params!(c::Cat{M, V}, p::V) where {M, V}
    c.p .= p
end

function Assignment(
        node_labels, edge_list::EdgeList{E},
        dist::Dist{D}) where {E, D <: Cat}
    dists = fit(dist, edge_list)
    realized = Matrix{Vector{Int}}(undef, n_groups, n_groups)
    counts = Matrix{Int}(undef, n_groups, n_groups)
    estimated = Matrix{Vector{Float64}}(undef, n_groups, n_groups)
    fill!(realized, zeros(Int, num_categories(unwrap(dist))))
    fill!(counts, 0)
    for u in 1:n_nodes
        g1 = node_labels[u]
        for (v, e) in iterate_neighbors(edge_list, u)
            g2 = node_labels[v]
            if u < v
                counts[minmax(g1, g2)...] += 1
                realized[minmax(g1, g2)...][e] += 1
            else
                break
            end
        end
    end
    for g2 in 1:n_groups
        for g1 in g2:n_groups
            estimated[g1, g2] = (counts[g1, g2] == 0) ?
                                zeros(Float64, num_categories(unwrap(dist))) :
                                (realized[g1, g2]) ./ counts[g1, g2]
        end
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
    w = WorkspaceDiscreteSwap{Dist{D},
        Matrix{Int}, Matrix{Float64}, Matrix{Float64},
        SymArray{D}, SymArray{Float64}}(
        deepcopy(θ), deepcopy(log_likelihood_per_group),
        counts, deepcopy(realized), deepcopy(estimated))
    return Assignment(
        node_labels, edge_list, dists, θ, log_likelihood_per_group, w)
end

function make_workspace(a::Assignment{E, Dist{D},
        F, W}) where {E, F, D <: Cat, W}
    return deepcopy(a.additional_workspace)
end

function make_swap!(ws::WorkspaceDiscreteSwap,
        a) where {E, F, D <: Categorical}
    ws.θ = deepcopy(a.θ)
    ws.log_likelihood_per_group = deepcopy(a.log_likelihood)
end
