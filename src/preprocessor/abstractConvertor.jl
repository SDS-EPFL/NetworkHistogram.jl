abstract type AbstractConvertor end

Base.broadcastable(o::T) where {T <: AbstractConvertor} = Ref(o)

"""
    Convert data from its original form to a processed form suitable for SBM estimation.
"""
function convert end

struct CategoricalConvertor{T} <: AbstractConvertor
    m::Int  # number of categories
    has_zero::Bool  # whether data contains zero values
    map::Dict{T, Int}
end

function num_bins(c::CategoricalConvertor)
    return c.m
end

function convert(c::CategoricalConvertor{T}, A::AbstractMatrix{T}) where {T}
    # Map original values to 1-based indices
    @error "to be implemented"
end

struct ContinuousConvertor{B, N, V <: AbstractVector{B}} <: AbstractConvertor
    zero_index::Int
    bins::V
end

function num_bins(c::ContinuousConvertor{B, N}) where {B, N}
    return N
end

## assume no singleton bins
function ContinuousConvertor(bins::AbstractVector{B}) where {B <:
                                                             Union{Interval, BareInterval}}
    bins = sort(bins, lt = lt = strictprecedes)
    N = length(bins) + 1
    zero_index = 1
    ContinuousConvertor{B, N, typeof(bins)}(zero_index, bins)
end

# assume bins are sorted and correctly cover the whole support
function (c::ContinuousConvertor{<:Union{Interval, BareInterval}})(x)
    iszero(x) && return c.zero_index
    x >= sup(c.bins[end]) && return length(c.bins) + 1
    x <= inf(c.bins[1]) && return c.zero_index + 1
    return findfirst(b -> in_interval(x, b), c.bins) + 1
end

function ContinuousConvertor(l, u, num_bins::Int)
    edges = collect(range(l, stop = u, length = num_bins + 1))
    bins = [bareinterval(edges[i], edges[i + 1]) for i in 1:num_bins]
    ContinuousConvertor(bins)
end
