"""
    struct ZeroInflatedCategorical{B, D} <: DiscreteUnivariateDistribution

A zero-inflated categorical distribution that combines a Bernoulli distribution with a categorical distribution.

# Fields
- `edge_proba::B`: The Bernoulli distribution representing the probability of zero.
- `dist::D`: The categorical distribution.

# Constructors
- `ZeroInflatedCategorical(p::Real, dist::D)`: Creates a zero-inflated categorical distribution with probability `p` of zero and categorical distribution `dist`.

# Mathematical Explanation
The zero-inflated categorical distribution modifies the original categorical distribution by introducing a probability `p` of zero. The `pmf` and `cdf` are adjusted accordingly:
- `pdf(x) = p * δ(x) + (1 - p) * pmf_original(x)`
- `cdf(x) = p * δ(x) + (1 - p) * cdf_original(x)`
where `δ(x)` is the Dirac delta function.
"""
struct ZeroInflatedCategorical{B, D} <: DiscreteUnivariateDistribution
    edge_proba::B
    dist::D
end

_dirac_delta(x) = x == 0 ? one(x) : zero(x)
_dirac_delta(x, lb, ub) = lb <= x <= ub ? one(x) : zero(x)

function ZeroInflatedCategorical(p::Real, dist::D) where {D}
    if p < 0
        p = zero(p)
    elseif p > 1
        p = one(p)
    end
    return ZeroInflatedCategorical(Bernoulli(1 - p), dist)
end

function ZeroInflatedCategorical(p::Real, probs::AbstractVector)
    if sum(probs) == 0
        probs_ = ones(length(probs)) / length(probs)
    else
        probs_ = probs / sum(probs)
    end
    if p < 0
        p = zero(p)
    elseif p > 1
        p = one(p)
    end
    return ZeroInflatedCategorical(p, Categorical(probs_))
end

function ZeroInflatedCategorical(vec_probs::AbstractVector)
    ZeroInflatedCategorical(vec_probs[1], vec_probs[2:end])
end

ZeroInflatedCategorical(k::Int) = ZeroInflatedCategorical(ones(k + 1) ./ (k + 1))

"""
    Distributions.pdf(d::ZeroInflatedCategorical, x::Real)

Computes the probability mass function (pmf) of the zero-inflated categorical distribution `d` at `x`.

# Mathematical Explanation
The `pmf` of the zero-inflated categorical distribution is given by:
- `pmf(x) = p * δ(x) + (1 - p) * pmf_original(x)`
where `p` is the probability of zero, `δ(x)` is the Dirac delta function, and `pmf_original(x)` is the pmf of the original categorical distribution.
"""
function Distributions.pdf(d::ZeroInflatedCategorical, x::Real)
    return pdf(d.edge_proba, zero(x)) * _dirac_delta(x) +
           pdf(d.edge_proba, one(x)) * pdf(d.dist, x)
end

"""
    rand(rng::Random.AbstractRNG, d::ZeroInflatedCategorical)

Generates a random sample from the zero-inflated categorical distribution `d` using the random number generator `rng`.
"""
function rand(rng::Random.AbstractRNG, d::ZeroInflatedCategorical)
    return rand(rng, d.edge_proba) * rand(rng, d.dist)
end

logpdf(d::ZeroInflatedCategorical, x::Real) = log(pdf(d, x))

minimum(d::ZeroInflatedCategorical) = min(minimum(d.dist), 0)

maximum(d::ZeroInflatedCategorical) = max(maximum(d.dist), 0)

insupport(d::ZeroInflatedCategorical, x::Real) = x == 0 || insupport(d.dist, x)

"""
    Distributions.cdf(d::ZeroInflatedCategorical, x::Real)

Computes the cumulative distribution function (cdf) of the zero-inflated categorical distribution `d` at `x`.

# Mathematical Explanation
The `cdf` of the zero-inflated categorical distribution is given by:
- `cdf(x) = p * δ(x) + (1 - p) * cdf_original(x)`
where `p` is the probability of zero, `δ(x)` is the Dirac delta function, and `cdf_original(x)` is the cdf of the original categorical distribution.
"""
function Distributions.cdf(d::ZeroInflatedCategorical, x::Real)
    return pdf(d.edge_proba, zero(x)) * _dirac_delta(x, 0, Inf) +
           pdf(d.edge_proba, one(x)) * cdf(d.dist, x)
end

function Distributions.params(d::ZeroInflatedCategorical)
    (first(params(d.edge_proba)), params(d.dist)...)
end

ncategories(d::ZeroInflatedCategorical) = ncategories(d.dist)

"""
    Distributions.fit(::Type{ZeroInflatedCategorical{B, D}}, data::AbstractArray, n_cat)

Fits a zero-inflated categorical distribution to the given data.
"""
function Distributions.fit(
        ::Type{ZeroInflatedCategorical{B, D}}, data::AbstractArray, n_cat) where {
        B, D <: Categorical}
    indices_0 = findall(x -> x == 0, data)
    p = length(indices_0) / length(data)
    if p != 1
        dist = fit_mle(Categorical, n_cat, data[setdiff(1:end, indices_0)])
        return ZeroInflatedCategorical(p, dist)
    else
        return ZeroInflatedCategorical(1.0, zeros(n_cat))
    end
end

function get_params_cat_like(dist::ZeroInflatedCategorical)
    p = first(params(dist.edge_proba))
    probs = vcat(params(dist.dist)...)
    return vcat(1 - p, probs .* p)
end

function Base.convert(::Type{<:ZeroInflatedCategorical}, d::D) where {D}
    return ZeroInflatedCategorical(1.0, d)
end

function Base.convert(T::Type{<:Categorical}, d::ZeroInflatedCategorical)
    return T(get_params_cat_like(d))
end
