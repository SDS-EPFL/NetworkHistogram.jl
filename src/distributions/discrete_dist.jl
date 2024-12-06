mutable struct DiscretizedDistribution{D, L} <: ContinuousUnivariateDistribution where {D, L}
    discretizer::D
    probs::L
end

function DiscretizedDistribution(d::D, n_bins::Int, support_bound = extrema(d)) where {D}
    disc = HybridDiscretizer(n_bins, support_bound..., 0.0)
    # for now we keep track of the non-edges as well
    probs = Distributions.Categorical(nlabels(discretizer))
    return DiscretizedDistribution(disc, probs)
end

function DiscretizedDistribution(discretizer::Discretizer)
    return DiscretizedDistribution(
        discretizer, Distributions.Categorical(nlabels(discretizer)))
end

function pdf(d::DiscretizedDistribution, x::Real)
    if !supports_encoding(d.discretizer, x)
        return 0.0
    end
    # for now suppose that the non-edges are encoded in the last bin
    bin = encode(d.discretizer, x)
    return pdf(d.probs, bin) / binwidth(d.discretizer)
end

function logpdf(d::DiscretizedDistribution, x::Real)
    if !supports_encoding(d.discretizer, x)
        return -Inf
    end
    # for now suppose that the non-edges are encoded in the last bin
    bin = encode(d.discretizer, x)
    return log(pdf(d.probs, bin)) - log(binwidth(d.discretizer))
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
    return supports_encoding(d.discretizer, x)
end

function Base.convert(::Type{DiscretizedDistribution}, d::D) where {D}
    return DiscretizedDistribution(d, 10)
end


function Distributions.ncategories(d::DiscretizedDistribution)
    return ncategories(d.probs)
end


function Distributions.fit(::Type{<:DiscretizedDistribution{D,L}},data) where {D,L}
    return fit(L, data)
end

function set_params!(d::DiscretizedDistribution{D, L}, params) where {D,L}
    d.probs = L(params)
end
