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

# Binary network
A = Symmetric(rand(0:1, 100, 100))
A[diagind(A)] .= 0

# Initial partition into 3 groups
initial_labels = rand(1:3, 100)

# Fit network histogram
params = GreedyParams()
result = nethist(A, Bernoulli(0.5), initial_labels, params)

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
