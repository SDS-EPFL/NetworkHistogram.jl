"""
    SBMEstimator

Abstract base type for all Stochastic Block Model (SBM) estimators.

All concrete estimator types should implement:
- `estimate(estimator, data, initial_labels; progress=true)`: Main estimation function
- `score(estimator)`: Return current objective value (if applicable)
"""
abstract type SBMEstimator end

"""
    SumGreedyEstimator{C, S, NodeR, StopR}

Greedy optimization estimator for Stochastic Block Models using sum-of-squares loss.

This estimator uses a greedy node-swapping algorithm to minimize the loss function:
    L = (1/n_edges) * Σᵢⱼ [count(i,j) - ||realized(i,j)||²/count(i,j)]

The algorithm iteratively swaps nodes between groups to improve the block model fit.

# Type Parameters
- `C`: Type for count matrices (usually symmetric array of integers)
- `S`: Type for realized value matrices (usually symmetric array of vectors)
- `NodeR <: NodeSwapRule`: Rule for selecting which nodes to swap
- `StopR <: StopRule`: Rule for determining when to stop optimization

# Fields
- `counts::C`: Number of possible edges between each pair of groups
- `counts_swap::C`: Working copy of counts for swap evaluation
- `realized::S`: Sum of observed edge values between each pair of groups
- `realized_swap::S`: Working copy of realized values for swap evaluation
- `max_iter::Int`: Maximum number of iterations
- `node_swap_rule::NodeR`: Strategy for selecting nodes to swap
- `stop_rule::StopR`: Criterion for early stopping

# Example
```julia
k = 5  # number of groups
counts = SymArray(k, 0)
counts_swap = SymArray(k, 0)
realized = SymArray(zero(SizedMatrix{k, k, MVector{m, Int}}))
realized_swap = SymArray(zero(SizedMatrix{k, k, MVector{m, Int}}))

estimator = SumGreedyEstimator(
    counts, counts_swap, realized, realized_swap,
    max_iter=100_000,
    node_swap_rule=RandomGroupSwap(),
    stop_rule=PreviousBestValue(1000, Inf, :min)
)

labels = estimate(estimator, data, initial_labels)
```
"""
struct SumGreedyEstimator{C, S, NodeR <: NodeSwapRule, StopR <: StopRule} <: SBMEstimator
    counts::C
    counts_swap::C
    realized::S
    realized_swap::S
    max_iter::Int
    node_swap_rule::NodeR
    stop_rule::StopR
end

"""
    score(estimator::SumGreedyEstimator)

Compute the current objective value (loss) for the estimator.

Lower values indicate better fit to a block model structure.
"""
function score(estimator::SumGreedyEstimator)
    return loss_function(estimator.realized, estimator.counts)
end

"""
    init!(estimator::SumGreedyEstimator, data, initial_labels)

Initialize the estimator's count and realized value matrices from data.

Iterates through the upper triangle of the adjacency matrix (i < j) to avoid
double-counting edges in undirected graphs. Updates both the main and swap
workspace matrices.

# Arguments
- `estimator::SumGreedyEstimator`: The estimator to initialize
- `data::AbstractMatrix`: Network adjacency matrix
- `initial_labels::Vector{Int}`: Initial group assignments for nodes
"""
function init!(estimator::SumGreedyEstimator, data, initial_labels)
    # Iterate over upper triangle to avoid double-counting edges
    @inbounds for j in axes(data, 2)
        label_j = initial_labels[j]
        for i in 1:(j - 1)  # More efficient than i < j check inside loop
            edge_value = data[i, j]
            if !isnothing(edge_value)
                label_i = initial_labels[i]

                # Update both main and swap workspaces
                add_realized(estimator.realized[label_i, label_j], edge_value)
                add_realized(estimator.realized_swap[label_i, label_j], edge_value)
                add_counts!(estimator.counts, edge_value, label_i, label_j)
                add_counts!(estimator.counts_swap, edge_value, label_i, label_j)
            end
        end
    end
end

"""
    estimate(estimator::SumGreedyEstimator, data, initial_labels; progress=true)

Estimate node group assignments using greedy optimization with node swapping.

# Algorithm
The algorithm proceeds as follows:
1. Initialize count and realized value matrices from data and initial labels
2. For each iteration:
   a. Select two nodes to swap according to the swap rule
   b. Tentatively swap them and update statistics
   c. Accept swap if it improves the loss, otherwise revert
   d. Check stopping criterion
3. Return final node labels

# Arguments
- `estimator::SumGreedyEstimator`: The estimator with configuration
- `data::AbstractMatrix`: Network adjacency matrix (n × n)
- `initial_labels::Vector{Int}`: Initial group assignments (length n)
- `progress::Bool`: Whether to show progress bar (default: true)

# Returns
- `node_labels::Vector{Int}`: Optimized group assignments for each node
"""
function estimate(estimator::SumGreedyEstimator, data, initial_labels; progress = true)
    # Initialize counts and realized values from data
    init!(estimator, data, initial_labels)
    initialise_stop_rule!(estimator.stop_rule, estimator)

    # Compute initial loss
    current_loss = score(estimator)

    # Start with initial labeling
    node_labels = copy(initial_labels)

    # Progress tracking
    pbar = ProgressUnknown(
        enabled = progress,
        showspeed = true,
        desc = "Greedy search: "
    )

    # Update progress bar only every N iterations to reduce overhead
    progress_update_interval = max(1, estimator.max_iter ÷ 1000)

    # Main optimization loop
    for iter in 1:(estimator.max_iter)
        # Select two nodes to potentially swap
        index1, index2 = select_indices_swap(node_labels, estimator.node_swap_rule)

        group1 = node_labels[index1]
        group2 = node_labels[index2]

        # Only process if nodes are in different groups
        if group1 != group2
            # Update swap workspace to reflect the proposed swap
            # Using @inbounds for performance - loop bounds are guaranteed safe
            @inbounds for j in axes(data, 1)
                # Skip the swapped nodes themselves
                if j == index1 || j == index2
                    continue
                end

                group_j = node_labels[j]
                edge_val_1 = data[j, index1]
                edge_val_2 = data[j, index2]

                # Update for node1: remove from group1, add to group2
                remove_realized(estimator.realized_swap[group1, group_j], edge_val_1)
                remove_counts!(estimator.counts_swap, edge_val_1, group1, group_j)
                add_realized(estimator.realized_swap[group2, group_j], edge_val_1)
                add_counts!(estimator.counts_swap, edge_val_1, group2, group_j)

                # Update for node2: remove from group2, add to group1
                remove_realized(estimator.realized_swap[group2, group_j], edge_val_2)
                remove_counts!(estimator.counts_swap, edge_val_2, group2, group_j)
                add_realized(estimator.realized_swap[group1, group_j], edge_val_2)
                add_counts!(estimator.counts_swap, edge_val_2, group1, group_j)
            end

            # Tentatively apply swap
            node_labels[index1] = group2
            node_labels[index2] = group1

            # Compute new loss
            new_loss = loss_function(estimator.realized_swap, estimator.counts_swap)

            # Accept or reject swap
            if new_loss < current_loss
                # Accept: commit swap to main workspace
                deepcopy!(estimator.realized, estimator.realized_swap)
                copy!(estimator.counts, estimator.counts_swap)
                current_loss = new_loss
            else
                # Reject: revert labels and workspace
                node_labels[index1] = group1
                node_labels[index2] = group2
                deepcopy!(estimator.realized_swap, estimator.realized)
                copy!(estimator.counts_swap, estimator.counts)
            end
        end

        # Update progress bar

        # Update progress bar only periodically to reduce overhead
        if progress && (iter % progress_update_interval == 0 || iter == estimator.max_iter)
            update!(
                pbar, iter;
                showvalues = [
                    ("loss", current_loss),
                    info_to_print(estimator.stop_rule)
                ])
        end

        # Check stopping criterion
        if stopping_rule(current_loss, estimator.stop_rule)
            @info "Stopping criterion met at iteration $iter"
            finish!(pbar)
            break
        end
    end

    return node_labels
end

"""
    loss_function(realized, counts)

Compute the normalized sum-of-squares loss for block model fitting.

The loss measures how well a block model fits the data by computing:
    L = (1/N) * Σᵢⱼ [count(i,j) - ||realized(i,j)||²/count(i,j)]

where the sum is over the upper triangle (i ≤ j) to avoid double-counting.

# Mathematical Interpretation
For each pair of groups (i,j):
- `count(i,j)` is the number of edges between groups i and j
- `realized(i,j)` is a vector of observed edge values
- The term `||realized(i,j)||²/count(i,j)` measures concentration of values
- Lower loss indicates better block structure (more homogeneous within blocks)

# Arguments
- `realized`: Symmetric array of realized edge value sums between groups
- `counts`: Symmetric array of edge counts between groups

# Returns
- Normalized loss value (lower is better)

# !warning
    This will need to be modified for other data types!
"""
@inline function loss_function(realized, counts::AbstractArray{<:Real})
    total_loss = 0.0
    total_edges = 0.0

    @inbounds for j in axes(realized, 2)
        for i in 1:j
            n_edges = counts[i, j]
            if n_edges > 0
                sum_squares = sum(abs2, realized[i, j])
                total_loss += n_edges - sum_squares / n_edges
                total_edges += n_edges
            end
        end
    end

    return total_edges > 0 ? total_loss / total_edges : 0.0
end

# this assumes that sum realized = counts
@inline function loss_function(realized, counts)
    total_loss = 0.0
    total_edges = 0.0
    @inbounds for j in axes(realized, 2)
        for i in 1:j
            for m in eachindex(realized[i, j])
                total_edges += realized[i, j][m]
                total_loss += realized[i, j][m] *
                              (1 -
                               _fast_div_(realized[i, j][m], counts[i, j][m]))
            end
        end
    end
    return total_loss / total_edges
end

@inline function _fast_div_(num::Real, denom::Real)
    num == 0.0 && denom == 0.0 && return 0.0
    return num / denom
end

# ============================================================================
# Count manipulation helpers
# ============================================================================

"""
    add_realized(parameter::AbstractArray, data_value::AbstractArray)

Add array data value to parameter array (for categorical edge values).
"""
@inline function add_realized(parameter::AbstractArray, data_value::AbstractArray)
    @inbounds parameter .+= data_value
end

"""
    remove_realized(parameter::AbstractArray, data_value::AbstractArray)

Remove array data value from parameter array (for categorical edge values).
"""
@inline function remove_realized(parameter::AbstractArray, data_value::AbstractArray)
    @inbounds parameter .-= data_value
end

"""
    add_realized(parameter::AbstractArray, data_value::Real)

Increment the count for a specific category (for categorical edge values).
"""
@inline function add_realized(parameter::AbstractArray, data_value::Real)
    @inbounds parameter[data_value] += 1
end

"""
    remove_realized(parameter::AbstractArray, data_value::Real)

Decrement the count for a specific category (for categorical edge values).
"""
@inline function remove_realized(parameter::AbstractArray, data_value::Real)
    @inbounds parameter[data_value] -= 1
end

@inline function add_counts!(
        counts::AbstractArray{T}, data_value::Real, group_i::Int, group_j::Int) where {T <:
                                                                                       Real}
    @inbounds counts[group_i, group_j] += one(T)
end

@inline function remove_counts!(
        counts::AbstractArray{T}, data_value::Real, group_i::Int, group_j::Int) where {T <:
                                                                                       Real}
    @inbounds counts[group_i, group_j] -= one(T)
end

@inline function add_counts!(
        counts::AbstractArray, data_value, group_i::Int, group_j::Int)
    @inbounds counts[group_i, group_j] .+= 1#data_value
end

@inline function remove_counts!(
        counts::AbstractArray, data_value, group_i::Int, group_j::Int)
    @inbounds counts[group_i, group_j] .-= 1#data_value
end

# ============================================================================
# Data preparation utilities
# ============================================================================

"""
    prepare_data_cat(A::AbstractMatrix{<:Real}, k; m=length(unique(A)), has_zero=zero(eltype(A)) in A)

Prepare categorical network data for SumGreedyEstimator.

Creates the necessary data structures (count matrices and realized value tensors)
for estimating a categorical Stochastic Block Model with k groups.

# Arguments
- `A::AbstractMatrix{<:Real}`: Adjacency matrix with categorical edge values
- `k::Int`: Number of groups to partition nodes into
- `m::Int`: Number of edge categories (default: inferred from unique values in A)
- `has_zero::Bool`: Whether the data contains zero values (default: auto-detected)

# Returns
A tuple containing:
- `data`: Preprocessed adjacency matrix (shifted if zero-indexed)
- `counts`: Symmetric k×k array for edge counts (initialized to 0)
- `counts_swap`: Workspace copy of counts for swap evaluation
- `realized`: Symmetric k×k array of m-dimensional count vectors (initialized to 0)
- `realized_swap`: Workspace copy of realized for swap evaluation

# Example
```julia
# Network with 3 edge types (0, 1, 2) for no edge, layer 1, layer 2
A = rand(0:2, 100, 100)
A = (A + A') .÷ 2  # Make symmetric

data, counts, counts_swap, realized, realized_swap = prepare_data_cat(A, k=5)
```

# Notes
- If data contains zeros, they are shifted to 1-indexing for categorical representation
- The realized arrays use StaticArrays.MVector for performance
- The symmetric array structure avoids redundant storage
"""
function prepare_data_cat(
        A::AbstractMatrix{<:Real},
        k::Int;
        m::Int = length(unique(A)),
        has_zero::Bool = zero(eltype(A)) in A
)
    @debug "Preparing data for categorical SBM with $m categories and $k groups."

    # Adjust data if zero-indexed (shift to 1-indexing for Julia)
    if has_zero
        @debug "Data contains zero values, using 1-based indexing."
        data = A .+ 1
    else
        data = A
    end

    # Initialize count matrices
    counts = SymArray(k, 0)
    counts_swap = SymArray(k, 0)

    # Initialize realized value tensors (k×k matrices of m-dimensional vectors)
    realized = SymArray(zero(SizedMatrix{k, k, MVector{m, Int}}))
    realized_swap = SymArray(zero(SizedMatrix{k, k, MVector{m, Int}}))

    return data, counts, counts_swap, realized, realized_swap
end
