"""
FastSymArray - Efficient symmetric matrix storage

This module provides `SymArray`, a memory-efficient storage for symmetric matrices
that only stores the upper triangle (including diagonal) of the matrix using a sparse matrix.
"""
module FastSymArray

using SparseArrays
using LinearAlgebra
import Base: eltype, convert, size, getindex, setindex!, copy!, similar,
             IndexStyle, axes, length, iterate, copyto!, fill!
import SparseArrays: getcolptr, nonzeros, FixedSparseCSC

export SymArray, eltype, deepcopy!, sum_tri_with_diag

"""
    SymArray{F} <: AbstractSparseMatrix{F, 2}

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

SymArray(::Type{F}, dims::Int...) where {F} = SymArray(F, dims)
function SymArray(::Type{F}, dims::NTuple{2, Int}) where {F}
    if dims[1] != dims[2]
        throw(ArgumentError("SymArray must be square, got dims=$(dims)"))
    end
    SymArray{F}(SparseMatrixCSC{F, Int}(make_csc_format(dims[1], F)...))
end

SymArray{F}(::UndefInitializer, dims::Int...) where {F} = SymArray{F}(undef, dims)
function SymArray{F}(::UndefInitializer, dims::NTuple{2, Int}) where {F}
    return SymArray(F, dims)
end

function make_csc_format(k::Int, ::Type{F}) where {F}
    k > 0 || throw(ArgumentError("Matrix dimension k=$k must be positive"))

    n_elements = div(k * (k + 1), 2)  # Number of non-zeros in upper triangle

    colptr = Vector{Int}(undef, k + 1)
    rowval = Vector{Int}(undef, n_elements)
    nzval = Vector{F}(undef, n_elements)

    @inbounds for j in 1:(k + 1)
        colptr[j] = div((j - 1) * j, 2) + 1
    end

    idx = 1
    @inbounds for j in 1:k
        for i in 1:j
            rowval[idx] = i
            idx += 1
        end
    end
    return k, k, colptr, rowval, nzval
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

# axes function
function axes(a::SymArray)
    return axes(a.uppertrian)
end

# length function
function length(a::SymArray)
    return length(a.uppertrian)
end

# faster indexing by avoiding search, modified from SparseArrays
Base.@propagate_inbounds function getindex(A::SymArray, i0::Integer, i1::Integer)
    i0, i1 = minmax(i0, i1)
    @boundscheck checkbounds(A, i0, i1)
    r1 = Int(@inbounds getcolptr(A.uppertrian)[i1])
    nonzeros(A.uppertrian)[r1 + i0 - 1]
end

# faster indexing by avoiding search, modified from SparseArrays
Base.@propagate_inbounds function setindex!(A::SymArray, v, i::Int, j::Int)
    i, j = minmax(i, j)
    @boundscheck checkbounds(A, i, j)
    r1 = Int(@inbounds getcolptr(A.uppertrian)[j])
    nonzeros(A.uppertrian)[r1 + i - 1] = v
end

function similar(a::SymArray, ::Type{T} = eltype(a), dims::Dims{2} = size(a)) where {T}
    return SymArray{T}(undef, dims)
end

function copy!(dest::SymArray{F}, src::SymArray{F}) where {F}
    size(dest) == size(src) || throw(DimensionMismatch("arrays must have the same size"))
    copy!(dest.uppertrian.nzval, src.uppertrian.nzval)
    return nothing
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
    return sum(a.uppertrian)
end

function convert(::Type{SymArray{F}}, a::AbstractMatrix{F}) where {F}
    @assert size(a, 1) == size(a, 2)
    k = size(a, 1)

    # Directly build upper triangle sparse matrix
    # Pre-allocate with exact size needed
    m, n, colptr, rowval, nzval = make_csc_format(k, F)
    idx = 1
    @inbounds for j in 1:k
        for i in 1:j
            nzval[idx] = a[i, j]
            idx += 1
        end
    end
    return SymArray(SparseMatrixCSC{F, Int}(m, n, colptr, rowval, nzval))
end

function deepcopy!(dest::SymArray{F}, src::SymArray{F}) where {F <: AbstractArray}
    dest_ = dest.uppertrian.nzval
    src_ = src.uppertrian.nzval
    @inbounds for index in eachindex(src_)
        if isassigned(dest_, index)
            copy!(dest_[index], src_[index])
        else
            dest_[index] = copy(src_[index])
        end
    end
    return dest
end

deepcopy!(dest::SymArray{F}, src::SymArray{F}) where {F <: Real} = copy!(dest, src)

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
    A = find_symarray(bc)
    return SymArray(similar(A.uppertrian, ElType))
end

# Custom similar for broadcasted SymArrays
function Base.similar(
        bc::Broadcast.Broadcasted{SymArrayStyle}, ::Type{Nothing})
    A = find_symarray(bc)
    return similar(Array{Nothing}, axes(bc))
end

# Helper function to find a SymArray in the broadcast tree
find_symarray(bc::Broadcast.Broadcasted) = find_symarray(bc.args)
find_symarray(args::Tuple) = find_symarray(args[1], Base.tail(args))
find_symarray(x) = x
find_symarray(args::Tuple{}) = nothing
find_symarray(a::SymArray, rest) = a
find_symarray(::Any, rest) = find_symarray(rest)

end
