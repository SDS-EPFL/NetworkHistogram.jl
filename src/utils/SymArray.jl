"""
FastSymArray - Efficient symmetric matrix storage

This module provides `SymArray`, a memory-efficient storage for symmetric matrices
that only stores the upper triangle (including diagonal) of the matrix using a sparse matrix.
"""
module FastSymArray

using SparseArrays
using LinearAlgebra
import Base: eltype, convert, size, getindex, setindex!, copy!, similar,
             IndexStyle, axes, length, iterate, copyto!
export SymArray, eltype, deepcopy!, sum_tri_with_diag

"""
    SymArray{F} <: AbstractArray{F, 2}

A symmetric matrix that stores only the upper triangle using a sparse matrix.

For a k×k symmetric matrix, only k(k+1)/2 elements are stored instead of k².
This implementation uses Julia's SparseMatrixCSC for efficient storage and access.

# Fields
- `uppertrian::SparseMatrixCSC{F, Int}`: Sparse matrix storing the upper triangle (i ≤ j)

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
mutable struct SymArray{F} <: AbstractSparseMatrix{F, Int}
    uppertrian::SparseMatrixCSC{F, Int}
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

    # Pre-allocate arrays with exact size needed for upper triangle
    n_elements = div(k * (k + 1), 2)
    I_indices = Vector{Int}(undef, n_elements)
    J_indices = Vector{Int}(undef, n_elements)
    values = Vector{F}(undef, n_elements)

    idx = 1
    for j in 1:k
        for i in 1:j
            I_indices[idx] = i
            J_indices[idx] = j
            values[idx] = deepcopy(d)
            idx += 1
        end
    end

    uppertrian = sparse(I_indices, J_indices, values, k, k)
    return SymArray{F}(uppertrian)
end

function SymArray(k::T, d::AbstractArray) where {T <: Real}
    k > 0 || throw(ArgumentError("Matrix dimension k=$k must be positive"))

    # Pre-allocate arrays with exact size needed for upper triangle
    n_elements = div(k * (k + 1), 2)
    I_indices = Vector{Int}(undef, n_elements)
    J_indices = Vector{Int}(undef, n_elements)
    values = Vector{typeof(d)}(undef, n_elements)

    idx = 1
    for j in 1:k
        for i in 1:j
            I_indices[idx] = i
            J_indices[idx] = j
            values[idx] = deepcopy(d)
            idx += 1
        end
    end

    uppertrian = sparse(I_indices, J_indices, values, k, k)
    return SymArray{typeof(d)}(uppertrian)
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

function size(a::SymArray)
    return size(a.uppertrian)
end

# IndexStyle trait - use CartesianIndex for 2D arrays
Base.IndexStyle(::Type{<:SymArray}) = IndexCartesian()

# axes function
function axes(a::SymArray)
    return axes(a.uppertrian)
end

# length function
function length(a::SymArray)
    return length(a.uppertrian)
end

Base.@propagate_inbounds function getindex(a::SymArray{F}, i::Int, j::Int) where {F}
    @boundscheck checkbounds(a, i, j)
    if i <= j
        @inbounds return a.uppertrian[i, j]
    else
        @inbounds return a.uppertrian[j, i]
    end
end

Base.@propagate_inbounds function setindex!(a::SymArray{F}, v, i::Int, j::Int) where {F}
    @boundscheck checkbounds(a, i, j)
    if i <= j
        @inbounds a.uppertrian[i, j] = v
    else
        @inbounds a.uppertrian[j, i] = v
    end
end

# similar function for creating similar arrays
function similar(a::SymArray{F}) where {F}
    k = size(a, 1)
    return SymArray(k, zero(F))
end

function similar(a::SymArray, ::Type{T}) where {T}
    k = size(a, 1)
    return SymArray(k, zero(T))
end

function similar(a::SymArray, ::Type{T}, dims::Dims{2}) where {T}
    dims[1] == dims[2] || throw(ArgumentError("SymArray must be square"))
    return SymArray(dims[1], zero(T))
end

function copyto!(dest::SymArray{F}, src::SymArray{F}) where {F}
    size(dest) == size(src) || throw(DimensionMismatch("arrays must have the same size"))
    copyto!(dest.uppertrian, src.uppertrian)
    return dest
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
    return sum(a.uppertrian.nzval)
end

function eltype(::SymArray{F}) where {F}
    return F
end

function convert(::Type{SymArray{F}}, a::AbstractMatrix{F}) where {F}
    @assert size(a, 1) == size(a, 2)
    k = size(a, 1)

    # Directly build upper triangle sparse matrix
    # Pre-allocate with exact size needed
    I_indices = Vector{Int}(undef, div(k * (k + 1), 2))
    J_indices = Vector{Int}(undef, div(k * (k + 1), 2))
    values = Vector{F}(undef, div(k * (k + 1), 2))

    idx = 1
    for j in 1:k
        for i in 1:j
            I_indices[idx] = i
            J_indices[idx] = j
            values[idx] = a[i, j]
            idx += 1
        end
    end

    uppertrian = sparse(I_indices, J_indices, values, k, k)
    return SymArray{F}(uppertrian)
end

function convert(::Type{AbstractMatrix{F}}, a::SymArray{F}) where {F}
    # Reconstruct full symmetric matrix from upper triangle
    # m = upper + upper' - Diagonal(upper) creates the full symmetric matrix
    m = a.uppertrian + transpose(a.uppertrian) -
        SparseArrays.spdiagm(0 => diag(a.uppertrian))
    return Matrix(m)
end

function copy!(dest::SymArray{F}, src::SymArray{F}) where {F <: Real}
    copyto!(dest, src)
    return dest
end

function deepcopy!(dest::SymArray{F}, src::SymArray{F}) where {F <: AbstractArray}
    @inbounds for index in eachindex(dest)
        copyto!(dest[index], src[index])
    end
    return dest
end

# Broadcasting support - custom style to maintain symmetric structure
struct SymArrayStyle <: Broadcast.AbstractArrayStyle{2} end
SymArrayStyle(::Val{2}) = SymArrayStyle()

Base.BroadcastStyle(::Type{<:SymArray}) = SymArrayStyle()

# When broadcasting with scalars or other styles, keep SymArrayStyle
Base.BroadcastStyle(::SymArrayStyle, ::Broadcast.DefaultArrayStyle{0}) = SymArrayStyle()
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::SymArrayStyle) = SymArrayStyle()

# When broadcasting with other arrays, use default array style
function Base.BroadcastStyle(::SymArrayStyle, ::Broadcast.DefaultArrayStyle)
    Broadcast.DefaultArrayStyle{2}()
end
function Base.BroadcastStyle(::Broadcast.DefaultArrayStyle, ::SymArrayStyle)
    Broadcast.DefaultArrayStyle{2}()
end

# When broadcasting between SymArrays, keep SymArrayStyle
Base.BroadcastStyle(::SymArrayStyle, ::SymArrayStyle) = SymArrayStyle()

# Custom similar for broadcasted SymArrays
function Base.similar(
        bc::Broadcast.Broadcasted{SymArrayStyle}, ::Type{ElType}) where {ElType}
    # For mutating functions that return Nothing, don't allocate a SymArray
    if ElType === Nothing
        # Find the first SymArray in the broadcast expression
        A = find_first_symarray(bc)
        # Return a similar array with the same element type as the input
        # This allows the broadcast to work but the result won't be used
        return similar(Array{ElType}, axes(bc))
    end
    # Find the first SymArray in the broadcast expression to get dimensions
    A = find_first_symarray(bc)
    return SymArray(size(A, 1), zero(ElType))
end

# Helper function to find a SymArray in the broadcast tree
find_first_symarray(bc::Broadcast.Broadcasted) = find_first_symarray(bc.args)
find_first_symarray(args::Tuple{}) = error("No SymArray found in broadcast")
find_first_symarray(args::Tuple) = find_first_symarray_in_args(args[1], Base.tail(args))

# Handle direct SymArray
find_first_symarray_in_args(x::SymArray, rest) = x
# Handle Extruded SymArray (from broadcasting)
find_first_symarray_in_args(x::Broadcast.Extruded{<:SymArray}, rest) = x.x
# Handle nested broadcasts
find_first_symarray_in_args(x::Broadcast.Broadcasted, rest) = find_first_symarray(x)
# Keep searching
find_first_symarray_in_args(x, rest) = find_first_symarray(rest)

# Custom copyto! for efficient broadcasting
function Base.copyto!(dest::SymArray, bc::Broadcast.Broadcasted{SymArrayStyle})
    # Broadcast only over the upper triangle for efficiency
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    bc′ = Broadcast.preprocess(dest, bc)

    # Only compute upper triangle
    k = size(dest, 1)
    @inbounds for j in 1:k
        for i in 1:j
            dest[i, j] = bc′[CartesianIndex(i, j)]
        end
    end
    return dest
end

# For broadcasting that returns Nothing (like with mutating functions)
function Base.copyto!(dest::AbstractArray, bc::Broadcast.Broadcasted{SymArrayStyle})
    # Fall back to default behavior
    Broadcast.materialize!(dest, bc)
end

@inline function throwdm(axdest, axsrc)
    throw(DimensionMismatch("destination axes $axdest are not compatible with source axes $axsrc"))
end

end
