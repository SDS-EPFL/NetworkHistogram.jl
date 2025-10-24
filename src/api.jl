"""
    nethist(data_input, dist_user, initial_node_labels, params::GreedyParams, zero_inflated::Bool = false)

Estimate a network histogram (stochastic block model) from network data.

This is the main entry point for fitting a network histogram to your data. It performs
preprocessing, optimization, and returns an Assignment representing the estimated model.

# Arguments
- `data_input`: Network data (adjacency matrix or EdgeList)
- `dist_user`: Reference distribution for edge values (e.g., Bernoulli, Categorical)
- `initial_node_labels`: Initial group assignment for nodes (vector of integers 1:k)
- `params::GreedyParams`: Optimization parameters
- `zero_inflated::Bool`: Whether to use zero-inflated version of distribution (default: false)

# Returns
- `Assignment`: The fitted network histogram with optimized node groups and parameters

# Throws
- `ArgumentError`: If input validation fails

# Examples
```julia
using NetworkHistogram, LinearAlgebra
import NetworkHistogram: nethist, GreedyParams

# Binary network
A = Symmetric(rand(0:1, 100, 100))
A[diagind(A)] .= 0

# Initial partition into 3 groups
initial_labels = rand(1:3, 100)

# Fit network histogram
result = nethist(A, Bernoulli(0.5), initial_labels,  GreedyParams())

# Extract results
block_matrix = result.θ
node_groups = result.node_labels
ll = loglikelihood(result)
```

See also: [`GreedyParams`](@ref), [`Assignment`](@ref), [`BlockModel`](@ref)
"""
function nethist(data_input, dist_user, initial_node_labels,
        params::GreedyParams, zero_inflated::Bool = false)
    # Input validation
    if data_input isa AbstractMatrix
        n_rows, n_cols = size(data_input)
        if n_rows != n_cols
            throw(ArgumentError("Adjacency matrix must be square, got size ($n_rows, $n_cols)"))
        end
        n = n_rows
    elseif data_input isa EdgeList
        n = number_nodes(data_input)
    else
        throw(ArgumentError("data_input must be an AbstractMatrix or EdgeList"))
    end

    if length(initial_node_labels) != n
        throw(ArgumentError("initial_node_labels length ($(length(initial_node_labels))) must match number of nodes ($n)"))
    end

    k = length(unique(initial_node_labels))
    if k < 1
        throw(ArgumentError("Must have at least one group, got $k groups"))
    end
    if k > n
        throw(ArgumentError("Number of groups ($k) cannot exceed number of nodes ($n)"))
    end

    if !all(x -> x isa Integer && 1 <= x <= k, initial_node_labels)
        throw(ArgumentError("initial_node_labels must contain integers in range 1:$k"))
    end

    if params.max_iter < 1
        throw(ArgumentError("max_iter must be positive, got $(params.max_iter)"))
    end

    return _nethist(
        data_input, dist_user, initial_node_labels, params, Val(zero_inflated))
end

# Internal implementation with compile-time zero-inflation flag
function _nethist(data_input, dist_user, initial_node_labels,
        params::GreedyParams, zero_inflated)
    @debug "preprocessing data"
    dist = get_ref_dist(dist_user, zero_inflated)
    g = preprocess_data(data_input, dist, zero_inflated)

    @debug "started optimization"
    out = greedy_optimize(g, initial_node_labels, params)

    @info "finished optimization with loglikelihood $(loglikelihood(out))"
    return postprocess(out)
end

# Helper functions for preprocessing

function get_ref_dist(dist::D, ::Val{true}) where {D}
    return Dist(ZeroInflated(dist))
end
function get_ref_dist(dist::D, ::Val{false}) where {D}
    return Dist(dist)
end

function preprocess_data(data, dist::Dist, zero_inflated)
    A = EdgeList(_fast_compressed_obs(dist, data, zero_inflated))
    return A, dist
end

function postprocess(out)
    return out
end

function nethist_discrete_edges(A, initial_node_labels, params::GreedyParams,
        k = length(unique(initial_node_labels)))
    data, counts_main, counts_swap, realized, realized_swap = prepare_data_cat(A, k)
    m = length(unique(data))
    es = SumGreedyEstimator(
        counts_main, counts_swap, realized, realized_swap,
        params.max_iter, params.swap_rule, params.stop_rule)
    node_labels = estimate(es, data, initial_node_labels)
    sizes = counts(node_labels) ./ length(node_labels)

    parameters = similar(es.realized)
    @inbounds for j in 1:k, i in 1:k
        parameters[i, j] = [es.realized[i, j][c] / es.counts[i, j] for c in 1:m]
    end
    model = DecoratedSBM(Categorical.(parameters), sizes)
    return NethistResult(node_labels, model)
end

function nethist_binary_edges(A, initial_node_labels, params::GreedyParams,
        k = length(unique(initial_node_labels)))
    data, counts_main, counts_swap, realized, realized_swap = prepare_data_cat(A, k)

    es = SumGreedyEstimator(
        counts_main, counts_swap, realized, realized_swap,
        params.max_iter, params.swap_rule, params.stop_rule)
    node_labels = estimate(es, data, initial_node_labels)
    sizes = counts(node_labels) ./ length(node_labels)

    θ = Matrix{Float64}(undef, k, k)
    @inbounds for j in 1:k, i in 1:k
        θ[i, j] = es.realized[i, j][2] / es.counts[i, j]
    end
    model = SBM(θ, sizes)

    return NethistResult(node_labels, model)
end

# functions for postprocessing

struct NethistResult{S}
    node_labels::Vector{Int}
    model::S
end

function NethistResult(a::Assignment)
    return NethistResult(copy(a.node_labels), to_block_model(a))
end

function to_block_model(a::Assignment{
        E, Dist{D}}) where {E, D <: Union{Bernoulli, Distributions.Bernoulli}}
    sizes = counts(a.node_labels) ./ length(a.node_labels)
    θ::Matrix{Float64} = map(x -> first(params(unwrap(x))), a.θ)
    return SBM(θ, sizes)
end

function to_block_model(a::Assignment)
    @info "Converting Assignment to DecoratedSBM"
    sizes = counts(a.node_labels) ./ length(a.node_labels)
    return DecoratedSBM(unwrap.(a.θ), sizes)
end

function node_labels_to_latents(node_labels::AbstractVector{Int}, sbm)
    return map(label -> _label_to_latent(label, sbm), node_labels)
end

function _label_to_latent(label::Int, sbm)
    return sbm.cumsize[label] - eps()
end

function align_res_true_latents!(res, a::Assignment, latents)
    perm = order_groups(a, latents)
    permute!(res.model, perm)
    res.node_labels .= map(x -> findfirst(==(x), perm), a.node_labels)
end

function permute!(sbm, perm)
    permuted_theta = copy(sbm.θ)
    sbm.θ .= permuted_theta[perm, perm]
    sbm.size .= sbm.size[perm]
    sbm.cumsize .= cumsum(sbm.size)
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
    return sort(
        1:k, by = x -> Tuple(get(counts[x], g, 0) for g in 1:k), rev = true)
end
