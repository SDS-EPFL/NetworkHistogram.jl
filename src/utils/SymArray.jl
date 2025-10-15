"""
FastSymArray - Efficient symmetric matrix storage

This module provides `SymArray`, a memory-efficient storage for symmetric matrices
that only stores the upper triangle (including diagonal) of the matrix.
"""
module FastSymArray

import Base: eltype, convert
export SymArray, eltype

"""
    SymArray{F} <: AbstractArray{F, 2}

A symmetric matrix that stores only the upper triangle to save memory.

For a k×k symmetric matrix, only k(k+1)/2 elements are stored instead of k².

# Fields
- `d::Dict{Tuple{Int, Int}, F}`: Dictionary storing (i,j) → value for i ≤ j
- `k::Int`: Dimension of the square matrix

# Examples
```julia
# Create a 3×3 symmetric matrix initialized with zeros
sym = SymArray(3, 0.0)

# Access elements (symmetric)
sym[1, 2] = 5.0
sym[2, 1]  # Returns 5.0

# Convert from regular matrix
A = [1 2 3; 2 4 5; 3 5 6]
sym = SymArray(A)
```

See also: [`sum_tri_with_diag`](@ref)
"""
mutable struct SymArray{F} <: AbstractArray{F, 2}
    d::Dict{Tuple{Int, Int}, F}
    k::Int
end

"""
    SymArray(k::Int, d::F)

Create a k×k symmetric matrix initialized with copies of value `d`.

# Arguments
- `k::Int`: Dimension of the matrix (must be positive)
- `d::F`: Initial value for all entries

# Example
```julia
sym = SymArray(5, 0.0)  # 5×5 matrix of zeros
```
"""
function SymArray(k::T, d::F) where {F, T <: Real}
    k > 0 || throw(ArgumentError("Matrix dimension k=$k must be positive"))
    return SymArray{F}(
        Dict{Tuple{Int, Int}, F}(minmax(i, j) => deepcopy(d) for i in 1:k
        for j in i:k),
        k)
end

function SymArray(k::T, d::AbstractArray) where {T <: Real}
    k > 0 || throw(ArgumentError("Matrix dimension k=$k must be positive"))
    return SymArray{typeof(d)}(
        Dict{Tuple{Int, Int}, typeof(d)}(minmax(i, j) => deepcopy(d)
        for i in 1:k
        for j in i:k),
        k)
end

"""
    SymArray(d::AbstractMatrix{F})

Create a SymArray from an existing matrix. The matrix should be symmetric.
Validates symmetry with a tolerance for floating-point errors.
"""
function SymArray(d::AbstractMatrix{F}) where {F}
    size(d, 1) == size(d, 2) || throw(ArgumentError(
        "Input matrix must be square, got size $(size(d))"))

    # Validate symmetry for floating-point types
    if F <: AbstractFloat
        k = size(d, 1)
        max_asymmetry = zero(F)
        for j in 1:k, i in 1:(j - 1)
            max_asymmetry = max(max_asymmetry, abs(d[i, j] - d[j, i]))
        end
        tol = sqrt(eps(F)) * maximum(abs, d)
        if max_asymmetry > tol
            @warn "Input matrix has asymmetry up to $max_asymmetry (tolerance: $tol). Using upper triangle."
        end
    end

    return convert(SymArray{F}, d)
end

function Base.size(a::SymArray)
    return (a.k, a.k)
end

Base.@propagate_inbounds function Base.getindex(a::SymArray, i, j)
    @boundscheck checkbounds(a, i, j)
    @inbounds return a.d[minmax(i, j)]
end

Base.@propagate_inbounds function Base.setindex!(a::SymArray, v, i, j)
    @boundscheck checkbounds(a, i, j)
    @inbounds a.d[minmax(i, j)] = v
end

"""
    sum_tri_with_diag(a::SymArray)

Efficiently sum all elements in the symmetric matrix (counting each off-diagonal once).

# Returns
- Sum of all unique elements in the symmetric matrix

# Note
This is more efficient than `sum(a)` because it only sums stored elements.
"""
function sum_tri_with_diag(a::SymArray)
    return sum(values(a.d))
end

function eltype(::SymArray{F}) where {F}
    return F
end

function convert(::Type{SymArray{F}}, a::AbstractMatrix{F}) where {F}
    @assert size(a, 1) == size(a, 2)
    k = size(a, 1)
    res = SymArray(k, a[1, 1])
    for j in axes(a, 2)
        for i in axes(a, 1)
            if i <= j
                res[i, j] = a[i, j]
            end
        end
    end
    return res
end

function convert(::Type{AbstractMatrix{F}}, a::SymArray{F}) where {F}
    k = a.k
    m = zeros(F, k, k)
    for i in 1:k
        for j in i:k
            m[i, j] = a[i, j]
        end
    end
    return m
end

end
