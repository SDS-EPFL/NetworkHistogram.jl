struct HistDistribution{B, P, P2, T} <: ContinuousUnivariateDistribution
    bins::B
    probs::P
    cum_probs::P2
    lower_bound::T
    upper_bound::T
end

Base.broadcastable(d::HistDistribution) = Ref(d)

params(d::HistDistribution) = (d.bins, d.probs)

eltype(::HistDistribution{F}) where {F} = F

function rand(rng::AbstractRNG, d::HistDistribution)
    u = rand(rng)
    bin_idx = searchsortedfirst(d.cum_probs, u)
    return uniform_in(rng, d.bins, bin_idx)
end

# assume bins are sorted and encoded by their upper bounds
function uniform_in(rng::AbstractRNG, bins::AbstractVector{B},
        bin_idx::Int) where {B <: Union{Interval, BareInterval}}
    l, u = inf(bins[bin_idx]), sup(bins[bin_idx])
    return rand(rng) * (u - l) + l
end

function HistDistribution(ps, c::ContinuousConvertor)
    cum_ps = cumsum(ps)
    intervals = vcat(bareinterval(0.0), c.bins)
    return HistDistribution(intervals, ps, cum_ps, inf(c.bins[1]), sup(c.bins[end]))
end

function HistDistribution(bins, ps)
    cum_ps = cumsum(ps)
    return HistDistribution(bins, ps, cum_ps, inf(bins[2]), sup(bins[end]))
end

### For Graphons compatibility
support(d::HistDistribution) = d.bins
_extract_param(d::HistDistribution, k) = d.probs[k]
# specialization for DiscreteNonParametric as it requires support to be specified
function convert_to_params(centers,
        sbm::DecoratedSBM{HistDistribution{B, P, P2, T}}) where {B, P, P2, T}
    s = support(sbm.θ[1, 1])
    return [HistDistribution(s, centers[:, i]) for i in axes(centers, 2)]
end

function logpdf(d::HistDistribution, x::Real)
    if x < d.lower_bound
        x = d.lower_bound + eps()
    elseif x >= d.upper_bound
        x = d.upper_bound - eps()
    end
    bin_idx = findfirst(b -> in_interval(x, b), d.bins) # potentially slow
    p = d.probs[bin_idx]
    # return sum(d.probs .^ 2) - p^2 + (1 - p)^2
    bin_idx == 1 && return log(p + eps())
    bin_width = sup(d.bins[bin_idx]) - inf(d.bins[bin_idx])
    return log(eps() + p / bin_width)
end
