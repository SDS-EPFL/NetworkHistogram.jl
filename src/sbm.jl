struct SBM{T, K, F <: Real} <: AbstractMatrix{T}
    sizes::Vector{F}
    probs::SymmetricTensor{T, K, 2}
end

function _check_sizes(sizes)
    @assert sum(sizes)≈1 "Sizes must sum to 1, got $(sum(sizes))"
    return sizes
end

function _check_sizes(sizes::Vector{Int})
    return sizes ./ sum(sizes)
end

function initialize_sbm(sizes::Vector, dist, k = length(sizes))
    sizes = _check_sizes(sizes)
    n_dims = binomial(k + 1, 2)
    probs = Vector{typeof(dist)}(undef, n_dims)
    fill!(probs, dist)
    return SBM(sizes, SymmetricTensor(probs, Val(k), Val(2)))
end

function initialize_sbm(sizes::GroupSize, dist, k = length(sizes))
    size_bins = sizes ./ sum(sizes)
    n_dims = binomial(k + 1, 2)
    probs = Vector{typeof(dist)}(undef, n_dims)
    fill!(probs, dist)
    return SBM(size_bins, SymmetricTensor(probs, Val(k), Val(2)))
end

function initialize_sbm(k::Int, dist)
    return initialize_sbm(ones(k) / k, dist)
end

number_blocks(::SBM{T, K}) where {T, K} = K

Base.size(s::SBM) = size(s.probs)
Base.ndims(::SBM) = 2
Base.eltype(::SBM{T, K}) where {T, K} = T
Base.setindex!(s::SBM, v, i, j) = setindex!(s.probs, v, i, j)
Base.@propagate_inbounds function Base.getindex(s::SBM, i, j)
    return getindex(s.probs, i, j)
end

function sample(
        rng::Random.AbstractRNG, sbm::SBM, node_labels::Vector{Int})
    n_nodes = length(node_labels)
    A = BitMatrix(undef, n_nodes, n_nodes)
    for i in 1:n_nodes
        A[i, i] = zero(eltype(A))
        for j in (i+1):n_nodes
            A[i, j] = Random.rand(rng, sbm[node_labels[i], node_labels[j]])
            A[j, i] = A[i, j]
        end
    end
    return sparse(A), node_labels
end

function sample(sbm::SBM, node_labels::Vector{Int})
    sample(Random.default_rng(), sbm, node_labels)
end
function sample(
        rng::Random.AbstractRNG, sbm::SBM, n_nodes::Int, sorted = true)
    n_blocks = number_blocks(sbm)
    node_labels = StatsBase.sample(
        rng, 1:n_blocks, StatsBase.weights(sbm.sizes), n_nodes, replace = true)
    if sorted
        sort!(node_labels)
    end
    return sample(rng, sbm, node_labels)
end

sample(sbm::SBM, n_nodes::Int) = sample(Random.default_rng(), sbm, n_nodes)
