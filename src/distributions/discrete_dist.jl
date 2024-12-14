mutable struct DiscretizedDistribution{D, L} <:
               ContinuousUnivariateDistribution where {D, L}
    discretizer::D
    probs::L
end

function DiscretizedDistribution(d::D, n_bins::Int, support_bound = extrema(d)) where {D}
    disc = DiscretizerZeroToZero(n_bins, support_bound...)
    ps = zeros(non_zero_labels_counts(disc))
    for i in 1:non_zero_labels_counts(disc)
        lb, ub = NetworkHistogram.decode(disc, i)
        ps[i] = cdf(d, ub) - cdf(d, lb)
    end
    probs = ZeroInflatedCategorical(0.0, ps)
    return DiscretizedDistribution(disc, probs)
end

function DiscretizedDistribution(d::ZeroInflated, n_bins::Int, support_bound = extrema(d))
    disc = DiscretizerZeroToZero(n_bins, support_bound...)
    ps = zeros(non_zero_labels_counts(disc))
    for i in 1:non_zero_labels_counts(disc)
        lb, ub = decode(disc, i)
        ps[i] = cdf(d, ub) - cdf(d, lb)
    end
    probs = ZeroInflatedCategorical(pdf(d, 0.0), ps)
    return DiscretizedDistribution(disc, probs)
end

function DiscretizedDistribution(discretizer::Discretizer)
    return DiscretizedDistribution(
        discretizer, ZeroInflatedCategorical(non_zero_labels_counts(discretizer)))
end

function rand(rng::Random.AbstractRNG, d::DiscretizedDistribution)
    bin = rand(rng, d.probs)
    return _decode_randomly(rng, d.discretizer, bin)
end

function minimum(d::DiscretizedDistribution)
    return minimum(d.discretizer)
end

function maximum(d::DiscretizedDistribution)
    return maximum(d.discretizer)
end

function insupport(d::DiscretizedDistribution, x::Real)
    return support_encoding(d.discretizer, x)
end

function Base.convert(::Type{DiscretizedDistribution}, d::D) where {D}
    return DiscretizedDistribution(d, 10)
end

function Distributions.ncategories(d::DiscretizedDistribution)
    return ncategories(d.probs)
end

function Distributions.fit(::Type{<:DiscretizedDistribution{D, L}}, data) where {D, L}
    return fit(L, data)
end

function set_params!(d::DiscretizedDistribution{D, L}, params) where {D, L}
    d.probs = L(params...)
end

# fast trick, will fail if discretizer put other categorical bins....
function pdf(d::DiscretizedDistribution, x::Real)
    if x == 0
        return pdf(d.probs, 0)
    end
    if !support_encoding(d.discretizer, x)
        return 0.0
    end
    bin = encode(d.discretizer, x)
    return pdf(d.probs, bin) / binwidth(d.discretizer)
end

function logpdf(d::DiscretizedDistribution, x::Real)
    if !support_encoding(d.discretizer, x)
        return -Inf
    end
    if x == 0
        return log(pdf(d.probs, 0))
    end
    bin = encode(d.discretizer, x)
    return log(pdf(d.probs, bin)) - log(binwidth(d.discretizer))
end

#lazy cdf computation, not efficient
function Distributions.cdf(
        d::DiscretizedDistribution{D, P}, x::Real) where {D, P <: ZeroInflatedCategorical}
    !insupport(d, x) && return 0.0
    bin = encode(d.discretizer, x)
    result = (x == 0) * cdf(d.probs, 0)
    if bin != 0
        lb, ub = decode(d.discretizer, bin)
        result += cdf(d.probs, bin - 1) +
                  (cdf(d.probs, bin) - cdf(d.probs, bin - 1)) * (x - lb) / (ub - lb)
    end
    return result
end
