"""
    struct ZeroInflated{B, D} <: ContinuousUnivariateDistribution

A zero-inflated distribution that combines a Bernoulli distribution with a continuous distribution.

# Fields
- `edge_proba::B`: The Bernoulli distribution representing the probability of zero.
- `dist::D`: The continuous distribution.

# Constructors
- `ZeroInflated(p::Real, dist::D)`: Creates a zero-inflated distribution with probability `p` of zero and continuous distribution `dist`.

# Mathematical Explanation
The zero-inflated distribution modifies the original distribution by introducing a probability `p` of zero. The `pdf` and `cdf` are adjusted accordingly:
- `pdf(x) = p * δ(x) + (1 - p) * pdf_original(x)`
- `cdf(x) = p * δ(x) + (1 - p) * cdf_original(x)`
where `δ(x)` is the Dirac delta function.
"""
struct ZeroInflated{B, D} <: ContinuousUnivariateDistribution
    edge_proba::B
    dist::D
end

function ZeroInflated(p::Real, dist::D) where {D}
    return ZeroInflated(Bernoulli(1 - p), dist)
end

"""
    Distributions.pdf(d::ZeroInflated, x::Real)

Computes the probability density function (pdf) of the zero-inflated distribution `d` at `x`.
"""
function Distributions.pdf(d::ZeroInflated, x::Real)
    return pdf(d.edge_proba, zero(x)) * _dirac_delta(x) +
           pdf(d.edge_proba, one(x)) * pdf(d.dist, x)
end

"""
    get_proba_zero(d::ZeroInflated)

Returns the probability of zero for the zero-inflated distribution `d`.
"""
function get_proba_zero(d::ZeroInflated)
    return pdf(d.edge_proba, 0)
end

"""
    rand(rng::Random.AbstractRNG, d::ZeroInflated)

Generates a random sample from the zero-inflated distribution `d` using the random number generator `rng`.
"""
function rand(rng::Random.AbstractRNG, d::ZeroInflated)
    return rand(rng, d.edge_proba) * rand(rng, d.dist)
end

logpdf(d::ZeroInflated, x::Real) = log(pdf(d, x))

minimum(d::ZeroInflated) = min(minimum(d.dist), 0)

maximum(d::ZeroInflated) = max(maximum(d.dist), 0)

insupport(d::ZeroInflated, x::Real) = x == 0 || insupport(d.dist, x)

"""
    Distributions.cdf(d::ZeroInflated, x::Real)

Computes the cumulative distribution function (cdf) of the zero-inflated distribution `d` at `x`.
"""
function Distributions.cdf(d::ZeroInflated, x::Real)
    return pdf(d.edge_proba, zero(x)) * _dirac_delta(x, zero(x), Inf) +
           cdf(d.dist, x) * pdf(d.edge_proba, one(x))
end

function Distributions.params(d::ZeroInflated)
    (first(params(d.edge_proba)), params(d.dist)...)
end

"""
    Distributions.fit(::Type{ZeroInflated{B, D}}, data::AbstractArray, n_cat)

Fits a zero-inflated distribution to the given data.
"""
function Distributions.fit(
        ::Type{ZeroInflated{B, D}}, data::AbstractArray, n_cat) where {B, D}
    indices_0 = findall(x -> x == 0, data)
    p = length(indices_0) / length(data)
    if p != 1
        return ZeroInflated(
            p, fit(D, data[setdiff(collect(eachindex(data)), indices_0)]))
    else
        return ZeroInflated(1.0, D())
    end
end
