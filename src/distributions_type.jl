struct Dist{D}
    dist::D
    counts::Int
    Dist(d,counts::Int) = counts < 0 ? error("Counts cannot be negative") : new{typeof(d)}(d, counts)
end



Dist(d) = Dist(d, 1)
zero(d::Dist) = Dist(zero(d.dist),0)

Base.broadcastable(x::Dist) = Ref(x)

function add_to(avgdist::Dist{D}, dist::D) where {D}
    return Dist(agg_params(avgdist.dist, dist, avgdist.counts / (avgdist.counts + 1), 1 / (avgdist.counts + 1)), avgdist.counts + 1)
end


function remove_from(avgdist::Dist{D}, dist::D) where {D}
    if avgdist.counts <= 0
        error("Cannot remove from a distribution with 0 counts")
    end
    # if avgdist.counts == 1 && params(avgdist) == params(dist)
    #     return Dist(zero(avgdist.dist), 0)
    # else
    #     error("Cannot remove from a distribution with 1 count unless the parameters are the same, got $(params(avgdist)) and $(params(dist))")
    # end
    return Dist(agg_params(avgdist.dist, dist, avgdist.counts / max(1,(avgdist.counts - 1)), - 1 / max(1,(avgdist.counts - 1))), avgdist.counts -1)
end


## probably this is fucked ...
add_to(d::Dist, dist::Dist) = add_to(d, dist.dist)
remove_from(d::Dist, dist::Dist) = remove_from(d, dist.dist)

for f in [:logpdf, :sample, :dist, :eltype, :params]
    @eval $f(d::Dist, args...) = $f(d.dist, args...)
end

fit(d::Dist, x) = Dist(fit(d.dist, x), d.counts)
loglikelihood(d::Dist, x) = sum(logpdf(d, y) for y in x)
unwrap(d::Dist) = d.dist


# expose compression step that assumes there is a pdf(d, typeof(compressed(x))) properly defined
# by default do nothing
_fast_compressed_obs(d, x) = x


# Bernoulli distribution (example)

struct Bernoulli{T<:Real}
    p::T
end


zero(d::Bernoulli) = Bernoulli(zero(d.p))
agg_params(d1::Bernoulli, d2::Bernoulli, w1, w2) = Bernoulli(w1 * d1.p + w2 * d2.p)
fit(::Bernoulli, x) = Bernoulli(mean(x))
dist(d1::Bernoulli, d2::Bernoulli) = abs(d1.p - d2.p)
logpdf(d::Bernoulli, x) = log(d.p * x + (1 - d.p) * (1 - x))
params(d::Bernoulli) = (d.p,)
eltype(d::Bernoulli) = Bool
sample(d::Bernoulli) = Bool(rand() <= d.p)
_fast_compressed_obs(d::Bernoulli, x) = x
