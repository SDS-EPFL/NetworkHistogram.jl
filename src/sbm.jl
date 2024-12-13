# TODO: remove BlockModel being a subtype of AbstractMatrix
# this was fun but useless and actually harmful

struct BlockModel{T, K, F <: Real} <: AbstractMatrix{T}
    sizes::Vector{F}
    probs::SymmetricTensor{T, K, 2}
end

function BlockModel(θ::AbstractMatrix{T}, sizes::Vector{F}) where {T, F <: Real}
    return BlockModel(sizes,
        SymmetricTensor([θ[i, j] for i in 1:size(θ, 1) for j in i:size(θ, 2)],
            Val(length(sizes)), Val(2)))
end

function edge_type(::BlockModel{T, K, F}) where {T, K, F}
    return eltype(T)
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
        rng::Random.AbstractRNG, sbm::BlockModel, node_labels::Vector{Int}, sorted = false)
    n_nodes = length(node_labels)
    if sorted
        sort!(node_labels)
    end
    A = zeros(edge_type(sbm), n_nodes, n_nodes)
    for j in 1:n_nodes
        for i in (j + 1):n_nodes
            A[i, j] = Random.rand(rng, sbm[node_labels[i], node_labels[j]])
        end
    end
    return sparse(Symmetric(A, :L)), node_labels
end


function draw_and_fill!(rng::Random.AbstractRNG, A, sbm::BlockModel, sorted = false)
    n_blocks = number_blocks(sbm)
    n_nodes = size(A, 1)
    node_labels = StatsBase.sample(
        rng, 1:n_blocks, StatsBase.weights(sbm.sizes), n_nodes, replace = true)
    if sorted
        sort!(node_labels)
    end
    @inbounds for j in 1:n_nodes
        for i in (j + 1):n_nodes
            A[i, j] = Random.rand(rng, sbm[node_labels[i], node_labels[j]])
        end
    end
    A .= Symmetric(A, :L)
end

draw_and_fill!(A, sbm, sorted = false) = draw_and_fill!(Random.default_rng(), A, sbm, sorted)

function sample(sbm::BlockModel, node_labels::Vector{Int}, sorted = false)
    sample(Random.default_rng(), sbm, node_labels, sorted)
end
function sample(
        rng::Random.AbstractRNG, sbm::BlockModel, n_nodes::Int, sorted = false)
    n_blocks = number_blocks(sbm)
    node_labels = StatsBase.sample(
        rng, 1:n_blocks, StatsBase.weights(sbm.sizes), n_nodes, replace = true)
    if sorted
        sort!(node_labels)
    end
    return sample(rng, sbm, node_labels)
end

function sample(sbm::BlockModel, n_nodes::Int, sorted = false)
    sample(Random.default_rng(), sbm, n_nodes, sorted)
end

function get_probability_matrix(sbm::BlockModel, node_labels::Vector{Int})
    return sbm.probs[node_labels, node_labels]
end

function _get_params_as_vec(dist::Distribution)
    return vcat(params(dist)...)
end


function latent_to_block_index(latents_vec, sbm::BlockModel)
    cum_sum_sizes = cumsum(sbm.sizes)
    cum_sum_sizes[end] = 1.0
    return [findfirst(x -> x >= l, cum_sum_sizes) for l in latents_vec]
end

"""
    best_alignment(fitted_sbm::BlockModel, true_sbm::BlockModel, tol = 0.01)

Find the best permutation of the blocks of `fitted_sbm` to match the blocks of `true_sbm` by
comparing the mean absolute difference of the parameters of the two models.
If the difference between the two models is less than `tol`, the function stops early.

!!! warning
    This function is not efficient for large numbers of blocks, as it uses brute force to
    find the best permutation.
"""
function best_alignment(fitted_sbm::BlockModel, true_sbm::BlockModel, tol = 0.01)
    k = number_blocks(fitted_sbm)
    if k != number_blocks(true_sbm)
        throw(ArgumentError("The number of blocks must be the same for both models"))
    end
    best_perm = nothing
    best_loss = Inf
    fitted_params = _get_params_as_vec.(fitted_sbm)
    true_params = _get_params_as_vec.(true_sbm)
    for perm in permutations(1:k)
        loss = sum(map(x -> sum(abs.(x)), fitted_params[perm] .- true_params))
        if loss < best_loss
            best_loss = loss
            best_perm = perm
        end
        if best_loss < tol
            break
        end
    end
    return best_perm
end

function align_sbm!(sbm::BlockModel, perm)
    sbm.probs .= sbm.probs[perm, perm]
    sbm.sizes .= sbm.sizes[perm]
end


"""
    order_groups(a::Assignment, latents::AbstractVector)

Order the groups of an assignment according to the true latents. This is an heuristic
approach, which is not guaranteed to find the true ordering of the groups.
"""
function order_groups(a::Assignment, latents::AbstractVector)
    n = number_nodes(a)
    k = number_groups(a)
    sort_perm = sortperm(latents)
    sorted_group_labels = a.node_labels[sort_perm]
    dummy_group_labels = repeat(1:k, inner = n ÷ k + 1)[1:n]
    counts = Dict(group => countmap(dummy_group_labels[sorted_group_labels .== group])
    for group in 1:k)
    return sort(1:k, by = x -> Tuple(get(counts[x], g, 0) for g in 1:k), rev = true)
end


function align_sbm_true_latents!(sbm::BlockModel, a::Assignment, latents)
    align_sbm!(sbm, order_groups(a, latents))
end
