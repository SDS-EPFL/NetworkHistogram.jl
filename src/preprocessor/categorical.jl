
# ============================================================================
# Data preparation utilities
# ============================================================================

"""
    prepare_data_cat(A::AbstractMatrix{<:Real}, k; m=length(unique(A)), has_zero=zero(eltype(A)) in A)

Prepare categorical network data for GreedyAverage.

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
