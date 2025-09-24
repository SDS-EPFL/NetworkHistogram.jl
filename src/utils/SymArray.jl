module FastSymArray

import Base: eltype, convert
export SymArray, eltype

mutable struct SymArray{F} <: AbstractArray{F, 2}
    d::Dict{Tuple{Int, Int}, F}
    k::Int
end

function SymArray(k::T, d::F) where {F, T <: Real}
    @assert k > 0
    return SymArray{F}(
        Dict{Tuple{Int, Int}, F}(minmax(i, j) => deepcopy(d) for i in 1:k
        for j in i:k),
        k)
end

function SymArray(k::T, d::AbstractArray) where {T <: Real}
    @assert k > 0
    return SymArray{typeof(d)}(
        Dict{Tuple{Int, Int}, typeof(d)}(minmax(i, j) => deepcopy(d)
        for i in 1:k
        for j in i:k),
        k)
end

function Base.size(a::SymArray)
    return (a.k, a.k)
end

Base.@propagate_inbounds function Base.getindex(a::SymArray, i, j)
    @boundscheck checkbounds(a, i, j)
    return a.d[minmax(i, j)]
end

Base.@propagate_inbounds function Base.setindex!(a::SymArray, v, i, j)
    @boundscheck checkbounds(a, i, j)
    a.d[minmax(i, j)] = v
end

function sum_tri_with_diag(a::SymArray)
    return sum(values(a.d))
end

function eltype(a::SymArray{F}) where {F}
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
