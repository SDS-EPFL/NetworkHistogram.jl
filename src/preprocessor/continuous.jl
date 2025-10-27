
### =======================================================================================
###                         [0,1] Continuous Convertor
### =======================================================================================

abstract type UnitIntervalConvertorType <: AbstractConvertor end

struct UnitIntervalConvertor{B <: AbstractVector} <: UnitIntervalConvertorType
    bins::B
end

function UnitIntervalConvertor(n::Int)
    zero_interval = Interval{:closed, :closed}(0.0, 0.0)
    edges = range(0.0, stop = 1.0, length = n + 1)
    bins = [Interval{:closed, :closed}(edges[i], edges[i + 1]) for i in 1:n]
    bins = vcat(zero_interval, bins)
    return UnitIntervalConvertor{typeof(bins)}(bins)
end

function num_bins(c::UnitIntervalConvertor)
    return length(c.bins)
end

function (c::UnitIntervalConvertor)(x::Real)
    return findfirst(b -> x ∈ b, c.bins)
end

function to_distribution(
        c::UnitIntervalConvertor, ps::AbstractVector{T}; kwargs...) where {T}
    @argcheck length(ps)==length(c.bins) "Length of probabilities must match number of bins"
    return HistDistribution(c.bins, SVector{length(ps), T}(ps))
end

# struct RegularUnitIntervalConvertor{N} <: UnitIntervalConvertorType
#     num_bins::Int
# end

### =======================================================================================
###                         Continuous Convertor
### =======================================================================================
# struct ContinuousConvertor{B, N, V <: AbstractVector{B}} <: AbstractConvertor
#     zero_index::Int
#     bins::V
# end

# function num_bins(c::ContinuousConvertor{B, N}) where {B, N}
#     return N
# end

# ## assume no singleton bins
# function ContinuousConvertor(bins::AbstractVector{B}) where {B <:
#                                                              Union{Interval, BareInterval}}
#     bins = sort(bins, lt = lt = strictprecedes)
#     N = length(bins) + 1
#     zero_index = 1
#     ContinuousConvertor{B, N, typeof(bins)}(zero_index, bins)
# end

# # assume bins are sorted and correctly cover the whole support
# function (c::ContinuousConvertor{<:Union{Interval, BareInterval}})(x)
#     iszero(x) && return c.zero_index
#     x >= sup(c.bins[end]) && return length(c.bins) + 1
#     x <= inf(c.bins[1]) && return c.zero_index + 1
#     return findfirst(b -> in_interval(x, b), c.bins) + 1
# end

# function ContinuousConvertor(l, u, num_bins::Int)
#     edges = collect(range(l, stop = u, length = num_bins + 1))
#     bins = [bareinterval(edges[i], edges[i + 1]) for i in 1:num_bins]
#     ContinuousConvertor(bins)
# end
