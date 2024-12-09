struct ZeroInflatedCategorical{B, D} <: DiscreteUnivariateDistribution
    edge_proba::B
    dist::D
end

_dirac_delta(x) = x == 0 ? one(x) : zero(x)

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

function Distributions.pdf(d::ZeroInflatedCategorical, x::Real)
    return pdf(d.edge_proba, 0) * _dirac_delta(x) + pdf(d.edge_proba, 1) * pdf(d.dist, x)
end

function rand(rng::Random.AbstractRNG, d::ZeroInflatedCategorical)
    return rand(rng, d.edge_proba) * rand(rng, d.dist)
end

logpdf(d::ZeroInflatedCategorical, x::Real) = log(pdf(d, x))

minimum(d::ZeroInflatedCategorical) = min(minimum(d.dist), 0)

maximum(d::ZeroInflatedCategorical) = max(maximum(d.dist), 0)

insupport(d::ZeroInflatedCategorical, x::Real) = x == 0 || insupport(d.dist, x)

function Distributions.cdf(d::ZeroInflatedCategorical, x::Real)
    return pdf(d.edge_proba, 0) * _dirac_delta(x) + pdf(d.edge_proba, 1) * cdf(d.dist, x)
end

function Distributions.params(d::ZeroInflatedCategorical)
    (first(params(d.edge_proba)), params(d.dist)...)
end

ncategories(d::ZeroInflatedCategorical) = ncategories(d.dist)

function Distributions.fit(
        ::Type{ZeroInflatedCategorical{B, D}}, data::AbstractArray, n_cat) where {B, D}
    indices_0 = findall(x -> x == 0, data)
    p = length(indices_0) / length(data)
    if p != 1
        dist = fit(D, data[setdiff(1:end, indices_0)])
        return ZeroInflatedCategorical(p, dist)
    else
        return ZeroInflatedCategorical(1.0, zeros(n_cat))
    end
end


function get_params_cat_like(dist::ZeroInflatedCategorical)
    p = first(params(dist.edge_proba))
    probs = vcat(params(dist.dist)...)
    return vcat(1-p, probs .* p)
end


function Base.convert(::Type{<:ZeroInflatedCategorical}, d::D) where {D}
    return ZeroInflatedCategorical(1.0, d)
end


function Base.convert(T::Type{<:Categorical}, d::ZeroInflatedCategorical)
    return T(get_params_cat_like(d))
end
