"""
    struct DiscretizedDistribution{D, L} <: ContinuousUnivariateDistribution

A discretized distribution that combines a discretizer with a zero-inflated categorical distribution.

# Fields
- `discretizer::D`: The discretizer used to discretize the continuous distribution.
- `probs::L`: The zero-inflated categorical distribution representing the discretized probabilities.

# Constructors
- `DiscretizedDistribution(d::D, n_bins::Int, support_bound = extrema(d))`: Creates a discretized distribution with `n_bins` bins and support bound `support_bound`.

# Mathematical Explanation
The discretized distribution modifies the original continuous distribution by dividing it into `n_bins` bins. The `pdf` and `cdf` are adjusted accordingly:
- `pdf(x) = pdf_discretized(bin) / bin_width`
- `cdf(x) = cdf_discretized(bin) + (cdf_discretized(bin + 1) - cdf_discretized(bin)) * progress_in_bin(x)`
"""
mutable struct DiscretizedDistribution{D, L} <:
               ContinuousUnivariateDistribution where {D, L}
    discretizer::D
    probs::L
end

function DiscretizedDistribution(
        d::D, n_bins::Int, support_bound = extrema(d)) where {D}
    disc = DiscretizerZeroToZero(n_bins, support_bound...)
    ps = zeros(non_zero_labels_counts(disc))
    for i in 1:non_zero_labels_counts(disc)
        lb, ub = decode(disc, i)
        ps[i] = cdf(d, ub) - cdf(d, lb)
    end
    probs = ZeroInflatedCategorical(0.0, ps)
    return DiscretizedDistribution(disc, probs)
end

function DiscretizedDistribution(
        d::ZeroInflated, n_bins::Int, support_bound = extrema(d))
    disc = DiscretizerZeroToZero(n_bins, support_bound...)
    ps = zeros(non_zero_labels_counts(disc))
    for i in 1:non_zero_labels_counts(disc)
        lb, ub = decode(disc, i)
        ps[i] = cdf(d, ub) - cdf(d, lb)
    end
    probs = ZeroInflatedCategorical(get_proba_zero(d), ps)
    return DiscretizedDistribution(disc, probs)
end

function DiscretizedDistribution(discretizer::Discretizer)
    return DiscretizedDistribution(
        discretizer, ZeroInflatedCategorical(non_zero_labels_counts(discretizer)))
end

"""
    rand(rng::Random.AbstractRNG, d::DiscretizedDistribution)

Generates a random sample from the discretized distribution `d` using the random number generator `rng`.
"""
function rand(rng::Random.AbstractRNG, d::DiscretizedDistribution)
    bin = rand(rng, d.probs)
    return _decode_randomly(rng, d.discretizer, bin)
end

minimum(d::DiscretizedDistribution) = minimum(d.discretizer)

maximum(d::DiscretizedDistribution) = maximum(d.discretizer)

function insupport(d::DiscretizedDistribution, x::Real)
    support_encoding(d.discretizer, x)
end

function Base.convert(::Type{DiscretizedDistribution}, d::D) where {D}
    return DiscretizedDistribution(d, 10)
end

ncategories(d::DiscretizedDistribution) = ncategories(d.probs)

function Distributions.fit(
        ::Type{<:DiscretizedDistribution{D, L}}, data) where {D, L}
    return fit(L, data)
end

function set_params!(d::DiscretizedDistribution{D, L}, params) where {D, L}
    d.probs = L(params...)
end

"""
    Distributions.pdf(d::DiscretizedDistribution, x::Real)

Computes the probability density function (pdf) of the discretized distribution `d` at `x`.

# Mathematical Explanation
The `pdf` of the discretized distribution is computed as:
- `pdf(x) = pdf_discretized(bin) / bin_width`
"""
function pdf(d::DiscretizedDistribution, x::Real)
    if x == 0
        return pdf(d.probs, zero(x))
    end
    if !support_encoding(d.discretizer, x)
        return zero(x)
    end
    bin = encode(d.discretizer, x)
    return pdf(d.probs, bin) / binwidth(d.discretizer)
end

"""
    Distributions.logpdf(d::DiscretizedDistribution, x::Real)

Computes the log of the probability density function (logpdf) of the discretized distribution `d` at `x`.
"""
function logpdf(d::DiscretizedDistribution, x::Real)
    if !support_encoding(d.discretizer, x)
        return -Inf
    end
    x == 0 && return log(pdf(d.probs, x))
    bin = encode(d.discretizer, x)
    return log(pdf(d.probs, bin)) - log(binwidth(d.discretizer))
end

"""
    Distributions.cdf(d::DiscretizedDistribution{D, P}, x::Real) where {D, P <: ZeroInflatedCategorical}

Computes the cumulative distribution function (cdf) of the discretized distribution `d` at `x`.

# Mathematical Explanation
The `cdf` of the discretized distribution is computed as:
- `cdf(x) = cdf_discretized(bin) + (cdf_discretized(bin + 1) - cdf_discretized(bin)) * progress_in_bin(x)`
"""
function Distributions.cdf(
        d::DiscretizedDistribution{D, P}, x::Real) where {
        D, P <: ZeroInflatedCategorical}
    x < minimum(d) && return zero(x)
    x > maximum(d) && return one(x)
    bin = encode(d.discretizer, x)
    result = (x == 0) * cdf(d.probs, x)
    if bin != 0
        result += cdf(d.probs, bin - 1) +
                  (cdf(d.probs, bin) - cdf(d.probs, bin - 1)) *
                  progress_in_bin(d.discretizer, x, bin)
    end
    return result
end
