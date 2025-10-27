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

export SymArray, eltype, deepcopy!, sum_tri_with_diag, make_sym_init

"""
    SymArray{F} <: AbstractSparseMatrix{F, 2}

A symmetric matrix that stores only the upper triangle using a sparse matrix.

For a k×k symmetric matrix, only k(k+1)/2 elements are stored instead of k².
This implementation uses Julia's SparseMatrixCSC for efficient storage and access.

# Fields
- `uppertrian::SparseMatrixCSC{F, Int}`: Sparse matrix storing the upper triangle (i ≤ j)

# Examples
```julia
# Create a 3×3 symmetric matrix
sym = SymArray{Float64}(undef, 3, 3)
sym .= 0.0

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

function make_sym_init(k, d::Real)
    a = SymArray{typeof(d)}(undef, k, k)
    fill!(a, d)
    return a
end

function make_sym_init(k, d)
    a = SymArray{typeof(d)}(undef, k, k)
    for j in 1:k, i in 1:j
        a[i, j] = deepcopy(d)
    end
    return a
end

@deprecate SymArray(k::Int, d::F) where {F} make_sym_init(k, d)

"""
    SymArray(d::AbstractMatrix{F})

Create a SymArray from an existing matrix.The matrix must be square and is assumed to be symmetric.
"""
function SymArray(d::AbstractMatrix{F}) where {F}
    m, n = size(d)
    m == n || throw(ArgumentError("Input matrix must be square, got size $(size(d))"))
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
    A.uppertrian.nzval[r1 + i0 - 1]
end

# faster indexing by avoiding search, modified from SparseArrays
Base.@propagate_inbounds function setindex!(A::SymArray, v, i::Int, j::Int)
    i, j = minmax(i, j)
    @boundscheck checkbounds(A, i, j)
    r1 = Int(@inbounds getcolptr(A.uppertrian)[j])
    A.uppertrian.nzval[r1 + i - 1] = v
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
    k, n = size(a)
    @assert k==n "Input matrix must be square, got size $(size(a))"

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
# struct SymArrayStyle <: Broadcast.AbstractArrayStyle{2} end
# SymArrayStyle(::Val{2}) = SymArrayStyle()

const SymArrayStyle = Broadcast.ArrayStyle{SymArray}

Base.BroadcastStyle(::Type{<:SymArray}) = Broadcast.ArrayStyle{SymArray}() # SymArrayStyle()

# When broadcasting with scalars or other styles, keep SymArrayStyle
Base.BroadcastStyle(::SymArrayStyle, ::Broadcast.DefaultArrayStyle{0}) = SymArrayStyle()
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::SymArrayStyle) = SymArrayStyle()

# When broadcasting with regular arrays (not scalars), defer to the array's style
# This ensures SymArray .+ Matrix returns Matrix, not SymArray
Base.BroadcastStyle(::SymArrayStyle, s::Broadcast.DefaultArrayStyle) = s
Base.BroadcastStyle(s::Broadcast.DefaultArrayStyle, ::SymArrayStyle) = s

# When broadcasting between SymArrays, keep SymArrayStyle
Base.BroadcastStyle(::SymArrayStyle, ::SymArrayStyle) = SymArrayStyle()

# Custom similar for broadcasted SymArrays
function Base.similar(
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SymArray}}, ::Type{ElType}) where {ElType}
    A = find_symarray(bc)
    if A == nothing
        return SymArray(similar(SparseMatrixCSC{ElType, Int}, axes(bc)...))
    else
        return SymArray(similar(A.uppertrian, ElType))
    end
end

# Helper function to find a SymArray in the broadcast tree
find_symarray(bc::Broadcast.Broadcasted) = find_symarray(bc.args)
find_symarray(args::Tuple) = find_symarray(args[1], Base.tail(args))
find_symarray(x) = x
find_symarray(args::Tuple{}) = nothing
find_symarray(a::SymArray, rest) = a
find_symarray(::Any, rest) = find_symarray(rest)

# Override broadcasted to eagerly evaluate when SymArrayStyle is involved
# This prevents issues with nested broadcasts losing the SymArray type
# hack, needs to be fixed later
function Broadcast.broadcasted(::SymArrayStyle, f, args...)
    # Eagerly materialize any nested Broadcasted{SymArrayStyle} to maintain type stability
    materialized_args = map(args) do arg
        if arg isa Broadcast.Broadcasted{SymArrayStyle}
            # Materialize nested SymArray broadcasts immediately
            return copy(arg)
        else
            return arg
        end
    end
    # Now create the broadcast with materialized args
    return Broadcast.Broadcasted{SymArrayStyle}(f, materialized_args)
end

# Specialized copyto! for efficient in-place broadcasting into SymArray
# This maintains the symmetric structure during broadcast operations
function Base.copyto!(dest::SymArray, bc::Broadcast.Broadcasted{SymArrayStyle})
    axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))

    _copyto_nzval!(dest, bc)
    return dest
    # # Try to use optimized nzval path for simple operations
    # if _can_use_nzval_broadcast(bc)
    #     return _copyto_nzval!(dest, bc)
    # end

    # # Fallback: iterate using CartesianIndices but only over upper triangle
    # bc′ = Broadcast.preprocess(dest, bc)
    # @inbounds for j in 1:size(dest, 2)
    #     for i in 1:j
    #         dest[i, j] = bc′[i, j]
    #     end
    # end
    # return dest
end

# Optimized copyto! that works directly on nzval arrays
function _copyto_nzval!(
        dest::SymArray{T}, bc::Broadcast.Broadcasted{SymArrayStyle}) where {T}
    # Replace SymArrays in the broadcast tree with their nzval arrays
    bc_nzval = _replace_with_nzval(bc)

    # Broadcast directly on the nzval array
    dest_nzval = nonzeros(dest.uppertrian)
    copyto!(dest_nzval, bc_nzval)

    return dest
end

# Replace SymArrays in broadcast tree with their nzval arrays
function _replace_with_nzval(bc::Broadcast.Broadcasted{SymArrayStyle})
    # Create new broadcasted with transformed arguments
    new_args = map(_replace_with_nzval, bc.args)
    # Don't specify style - let it be inferred
    return Broadcast.Broadcasted(bc.f, new_args)
end

function _replace_with_nzval(sa::SymArray)
    return sa.uppertrian.nzval
end

function _replace_with_nzval(bc::Broadcast.Broadcasted)
    # Recursively process nested broadcasts
    new_args = map(_replace_with_nzval, bc.args)
    return Broadcast.Broadcasted(bc.f, new_args)
end

function _replace_with_nzval(x)
    # For scalars and other types, return as-is
    return x
end

end
