"""
    BlockModel{D, V, M <: AbstractMatrix{D}}

A stochastic block model representation for network generation and analysis.

A block model is a piecewise constant graphon approximation where nodes are divided
into K blocks, and edges between blocks follow specific distributions.

# Fields
- `_dists::M`: Symmetric K×K matrix of edge distributions between blocks
- `sizes::V`: Proportions of nodes in each block (sums to 1.0)
- `cum_sizes::V`: Cumulative proportions for mapping latent variables to blocks

# Type Parameters
- `D`: Distribution type for edges (e.g., Bernoulli, Categorical, etc.)
- `V`: Vector type for storing proportions
- `M`: Matrix type for storing distributions

# Constructors

```julia
# Uniform block sizes with k blocks
BlockModel(k::Int, d::D)

# From an Assignment
BlockModel(a::Assignment)

# From node labels and parameter matrix
BlockModel(node_labels, θ)

# From a distribution matrix (infers uniform block sizes)
BlockModel(θ::AbstractMatrix)
```

# Examples
```julia
# Create a 3-block model with Bernoulli edges
bm = BlockModel(3, Bernoulli(0.5))

# Sample a network from the block model
latents, A = sample(bm, 100)  # 100 nodes

# Access block-to-block distribution
dist_12 = bm[1, 2]

# Map a latent variable to a block
block = map_ξ_to_block(bm, 0.3)
```

See also: [`Assignment`](@ref), [`sample`](@ref), [`get_probability_matrix`](@ref)
"""
struct BlockModel{D, V, M <: AbstractMatrix{D}}
    _dists::M
    sizes::V
    cum_sizes::V
end

"""
    BlockModel(k::Int, d::D) where {D}

Create a block model with `k` uniform-sized blocks, each initialized with distribution `d`.
"""
function BlockModel(k::Int, d::D) where {D}
    sizes = fill(1 / k, k)
    cumulative_sizes = cumsum(sizes)
    _dists = SymArray(k, d)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

"""
    BlockModel(a::Assignment)

Create a BlockModel from an Assignment, extracting the block proportions and
fitted distributions.
"""
function BlockModel(a::Assignment)
    k = length(unique(a.node_labels))
    sizes = proportions(a)
    cumulative_sizes = cumsum(sizes)
    _dists = unwrap.(a.θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

"""
    BlockModel(nodes_labels, θ)

Create a BlockModel from node labels and a distribution matrix θ.
"""
function BlockModel(nodes_labels, θ)
    k = length(unique(nodes_labels))
    sizes = counts(nodes_labels) / length(nodes_labels)
    cumulative_sizes = cumsum(sizes)
    _dists = unwrap.(θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

"""
    BlockModel(θ::AbstractMatrix{D}) where {D}

Create a BlockModel from a distribution matrix, assuming uniform block sizes.
"""
function BlockModel(θ::AbstractMatrix{D}) where {D}
    k = size(θ, 1)
    sizes = fill(1 / k, k)
    cumulative_sizes = cumsum(sizes)
    _dists = convert(SymArray{D}, θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

"""
    map_ξ_to_block(bm::BlockModel, ξ::Real)

Map a latent variable ξ ∈ [0,1] to its corresponding block index.

# Arguments
- `bm::BlockModel`: The block model
- `ξ::Real`: Latent variable in [0, 1]

# Returns
- `Int`: Block index (1 to k)
"""
function map_ξ_to_block(bm::BlockModel, ξ::T) where {T <: Real}
    return findfirst(x -> x >= ξ, bm.cum_sizes)
end

"""
    sample(bm::BlockModel, latents::Int, args...)

Sample a network from the block model by first generating `latents` random latent
variables, then sampling edges according to the block distributions.

# Arguments
- `bm::BlockModel`: The block model to sample from
- `latents::Int`: Number of nodes to generate
- `args...`: Additional arguments passed to edge sampling

# Returns
- Tuple of (latent_assignments, adjacency_matrix)
"""
function sample(bm::BlockModel, latents::Int, args...)
    latents = map(x -> map_ξ_to_block(bm, x), rand(latents))
    return latents, sample(bm, latents, args...)
end

"""
    sample(bm::BlockModel, latents::Vector, args...)

Sample a network from the block model given specific latent block assignments.

# Arguments
- `bm::BlockModel`: The block model to sample from
- `latents::Vector`: Block assignments for each node
- `args...`: Additional arguments passed to edge sampling

# Returns
- Adjacency matrix with sampled edges
"""
function sample(bm::BlockModel, latents::Vector{T}, args...) where {T}
    A = Array{eltype(bm[1, 1]), 2}(undef, length(latents), length(latents))
    for j in 1:length(latents)
        for i in 1:(j - 1)
            A[i, j] = A[j, i]
        end
        for i in (j + 1):length(latents)
            A[i, j] = sample(bm[latents[i], latents[j]], args...)
            A[j, i] = A[i, j]
        end
    end
    # Fill diagonal with zeros (no self-loops)
    for i in 1:length(latents)
        A[i, i] = zero(A[1, 2])
    end
    return A
end

# Base interface implementations for BlockModel

function Base.getindex(s::BlockModel, i::Int, j::Int)
    return s._dists[i, j]
end

function Base.setindex!(s::BlockModel, v, i::Int, j::Int)
    s._dists[i, j] = v
end

function Base.size(s::BlockModel)
    return (s._dists.k, s._dists.k)
end

"""
    getindex(bm::BlockModel, i::Real, j::Real)

Index into the block model using latent variables ξᵢ, ξⱼ ∈ [0,1].

Maps latent variables to their corresponding blocks and returns the
distribution between those blocks.
"""
function Base.getindex(s::BlockModel, i::Real, j::Real)
    k = findfirst(x -> x ≥ i, s.cum_sizes)
    l = findfirst(x -> x ≥ j, s.cum_sizes)
    return s._dists[k, l]
end

function Base.setindex!(s::BlockModel, v, i::Real, j::Real)
    k = findfirst(x -> x ≥ i, s.cum_sizes)
    l = findfirst(x -> x ≥ j, s.cum_sizes)
    s._dists[k, l] = v
end

"""
    ordered_latents(bm::BlockModel, n::Int)

Generate `n` ordered (sorted) latent block assignments from the block model.

# Returns
- Sorted vector of block assignments
"""
function ordered_latents(bm::BlockModel, n::Int)
    return sort(map(x -> map_ξ_to_block(bm, x), rand(n)))
end

"""
    get_probability_matrix(bm::BlockModel, latents::AbstractVector, default_dist=nothing)

Generate a node-level probability matrix from a block model and latent assignments.

Creates an n×n matrix where entry (i,j) contains the distribution for the edge
between nodes i and j, based on their block assignments.

# Arguments
- `bm::BlockModel`: The block model
- `latents::AbstractVector`: Block assignment for each node
- `default_dist`: Distribution for diagonal entries (defaults to zero(bm[1,1]) if not provided)

# Returns
- `Matrix`: n×n matrix of distributions

# Example
```julia
bm = BlockModel(3, Bernoulli(0.5))
latents = [1, 1, 2, 2, 3]
prob_matrix = get_probability_matrix(bm, latents)
```
"""
function get_probability_matrix(
        bm::BlockModel{D}, latents::AbstractVector, default_dist = nothing) where {D}
    # Set default distribution for diagonal (no self-loops)
    if isnothing(default_dist)
        try
            default_dist = zero(bm[1, 1])
        catch e
            if !isa(e, MethodError)
                rethrow(e)
            end
            error("Please provide a default distribution for the diagonal as it could not be inferred")
        end
    end
    n = length(latents)
    A = Array{D, 2}(undef, n, n)
    for j in 1:n
        for i in 1:n
            if i == j
                A[i, i] = default_dist
            else
                A[i, j] = bm[latents[i], latents[j]]
            end
        end
    end
    return A
end

"""
    get_probability_matrix(a::Assignment, default_dist=nothing, node_labels=a.node_labels)

Generate a node-level probability matrix from an Assignment.

# Arguments
- `a::Assignment`: The assignment
- `default_dist`: Distribution for diagonal entries (default: nothing)
- `node_labels`: Custom node labels to use (default: a.node_labels)

# Returns
- `Matrix`: Probability matrix based on the assignment's block structure
"""
function get_probability_matrix(
        a::Assignment, default_dist = nothing, node_labels = a.node_labels)
    return get_probability_matrix(BlockModel(a.θ), node_labels, default_dist)
end

"""
    align_sbm!(sbm::BlockModel, perm)

Permute the blocks of a stochastic block model according to permutation `perm`.

This modifies the block model in-place, reordering blocks and updating the
cumulative sizes accordingly.

# Arguments
- `sbm::BlockModel`: The block model to modify (modified in-place)
- `perm`: Permutation vector for reordering blocks
"""
function align_sbm!(sbm::BlockModel, perm)
    sbm._dists .= sbm._dists[perm, perm]
    sbm.sizes .= sbm.sizes[perm]
    sbm.cum_sizes .= cumsum(sbm.sizes)
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

function align_sbm_true_latents!(sbm::BlockModel, a::Assignment, latents)
    align_sbm!(sbm, order_groups(a, latents))
end
