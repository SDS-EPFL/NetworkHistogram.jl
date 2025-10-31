struct HistDistribution{B, P, P2} <: ContinuousUnivariateDistribution
    bins::B
    probs::P
    cum_probs::P2
end

Base.broadcastable(d::HistDistribution) = Ref(d)

params(d::HistDistribution) = (d.bins, d.probs)

function rand(rng::AbstractRNG, d::HistDistribution)
    u = rand(rng)
    bin_idx = searchsortedfirst(d.cum_probs, u)
    return rand(rng, d.bins[bin_idx])
end

function HistDistribution(bins, ps)
    cum_ps = SVector(cumsum(ps)...)
    return HistDistribution{typeof(bins), typeof(ps), typeof(cum_ps)}(bins, ps, cum_ps)
end

function logpdf(d::HistDistribution, x::Real)
    # potentially slow
    bin_idx = findfirst(b -> x ∈ b, d.bins)
    p = d.probs[bin_idx]
    bin_idx == 1 && return log(p)
    return log(p) - log(width(d.bins[bin_idx]))
end

function pdf(d::HistDistribution, x::Real)
    # potentially slow
    bin_idx = findfirst(b -> x ∈ b, d.bins)
    p = d.probs[bin_idx]
    bin_idx == 1 && return p
    return p / width(d.bins[bin_idx])
end

### For Graphons compatibility
support(d::HistDistribution) = d.bins
_extract_param(d::HistDistribution, k) = d.probs[k]

function convert_to_params(centers,
        sbm::DecoratedSBM{HistDistribution{B, P, P2}}) where {B, P, P2}
    s = sbm.θ[1, 1].bins
    return [HistDistribution(s, convert(P, centers[:, i])) for i in axes(centers, 2)]
end
