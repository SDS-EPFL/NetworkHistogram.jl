struct ZeroInflated{B, D} <: ContinuousUnivariateDistribution
    edge_proba::B
    dist::D
end

function ZeroInflated(p::Real, dist::D) where {D}
    return ZeroInflated(Bernoulli(1 - p), dist)
end

function Distributions.pdf(d::ZeroInflated, x::Real)
    return pdf(d.edge_proba, 0) * _dirac_delta(x) + pdf(d.edge_proba, 1) * pdf(d.dist, x)
end

function get_proba_zero(d::ZeroInflated)
    return pdf(d.edge_proba, 0)
end

function rand(rng::Random.AbstractRNG, d::ZeroInflated)
    return rand(rng, d.edge_proba) * rand(rng, d.dist)
end

logpdf(d::ZeroInflated, x::Real) = log(pdf(d, x))

minimum(d::ZeroInflated) = min(minimum(d.dist), 0)

maximum(d::ZeroInflated) = max(maximum(d.dist), 0)

insupport(d::ZeroInflated, x::Real) = x == 0 || insupport(d.dist, x)

function Distributions.cdf(d::ZeroInflated, x::Real)
    return pdf(d.edge_proba, 0) * _dirac_delta(x, 0, Inf) +
           cdf(d.dist, x) * pdf(d.edge_proba, 1)
end

function Distributions.params(d::ZeroInflated)
    (first(params(d.edge_proba)), params(d.dist)...)
end

function Distributions.fit(
        ::Type{ZeroInflated{B, D}}, data::AbstractArray, n_cat) where {B, D}
    indices_0 = findall(x -> x == 0, data)
    p = length(indices_0) / length(data)
    if p != 1
        return ZeroInflated(p, fit(D, data[setdiff(collect(eachindex(data)), indices_0)]))
    else
        return ZeroInflated(1.0, D())
    end
end
