struct BlockModel{T, K, F <: Real} <: AbstractMatrix{T}
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
    return BlockModel(sizes, SymmetricTensor(probs, Val(k), Val(2)))
end

function initialize_sbm(sizes::GroupSize, dist, k = length(sizes))
    size_bins = sizes ./ sum(sizes)
    n_dims = binomial(k + 1, 2)
    probs = Vector{typeof(dist)}(undef, n_dims)
    fill!(probs, dist)
    return BlockModel(size_bins, SymmetricTensor(probs, Val(k), Val(2)))
end

function initialize_sbm(k::Int, dist)
    return initialize_sbm(ones(k) / k, dist)
end

number_blocks(::BlockModel{T, K, F}) where {T, K, F} = K

Base.size(s::BlockModel) = size(s.probs)
Base.ndims(::BlockModel) = 2
Base.eltype(::BlockModel{T, K, F}) where {T, K, F} = T
Base.setindex!(s::BlockModel, v, i, j) = setindex!(s.probs, v, i, j)
Base.@propagate_inbounds function Base.getindex(s::BlockModel, i, j)
    return getindex(s.probs, i, j)
end

function sample(
        rng::Random.AbstractRNG, sbm::BlockModel, node_labels::Vector{Int})
    n_nodes = length(node_labels)
    type_input = eltype(sbm.probs[1, 1])
    A = Matrix{type_input}(undef, n_nodes, n_nodes)
    for i in 1:n_nodes
        A[i, i] = zero(eltype(A))
        for j in (i + 1):n_nodes
            A[i, j] = Random.rand(rng, sbm[node_labels[i], node_labels[j]])
            A[j, i] = A[i, j]
        end
    end
    return sparse(A), node_labels
end

function sample(sbm::BlockModel, node_labels::Vector{Int})
    sample(Random.default_rng(), sbm, node_labels)
end
function sample(
        rng::Random.AbstractRNG, sbm::BlockModel, n_nodes::Int, sorted = true)
    n_blocks = number_blocks(sbm)
    node_labels = StatsBase.sample(
        rng, 1:n_blocks, StatsBase.weights(sbm.sizes), n_nodes, replace = true)
    if sorted
        sort!(node_labels)
    end
    return sample(rng, sbm, node_labels)
end

function sample(sbm::BlockModel, n_nodes::Int)
    sample(Random.default_rng(), sbm, n_nodes)
end
