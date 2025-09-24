
struct ZeroInflated{D, F}
    dist::D
    proba_zero::F
end

struct SampleZI{F}
    value::F
    iszero::Bool
end

function ZeroInflated(dist)
    return ZeroInflated(dist, 0.0)
end

function logpdf(zi::ZeroInflated{D, F}, x::SampleZI) where {D, F}
    if x.iszero
        return log(zi.proba_zero)
    else
        return log(1 - zi.proba_zero) + logpdf(zi.dist, x.value)
    end
end

# function logpdf(zi::ZeroInflated{D, F}, x) where {D, F}
#     if iszero(x)
#         return log(zi.proba_zero)
#     else
#         return log(1 - zi.proba_zero) + logpdf(zi.dist, x.value)
#     end
# end

function agg_params(
        zi1::ZeroInflated{D, F}, zi2::ZeroInflated{D, F}, w1, w2) where {D, F}
    new_proba_zero = w1 * zi1.proba_zero + w2 * zi2.proba_zero
    return ZeroInflated(
        agg_params(zi1.dist, zi2.dist, w1, w2),
        new_proba_zero)
end

zero(zi::ZeroInflated) = ZeroInflated(zero(zi.dist), 0.0)

eltype(zi::ZeroInflated{D, F}) where {D, F} = SampleZI{eltype(D)}
params(zi::ZeroInflated{D, F}) where {D, F} = (params(zi.dist)..., zi.proba_zero)

function fit(zi::ZeroInflated{D, F}, x::SampleZI) where {D, F}
    if x.iszero
        return ZeroInflated(zero(zi.dist), 1.0)
    else
        return ZeroInflated(fit(zi.dist, x.value), 0.0)
    end
end

function _fast_compressed_obs(zi::ZeroInflated, x, filter = iszero)
    return SampleZI(_fast_compressed_obs(zi.dist, x), filter(x))
end

function unwrap(d::Dist{ZeroInflated{B, D}}) where {B, D}
    #yeah I know again...
    return d.dist.dist
end

function get_proportion_observed(d::Dist{ZeroInflated{B, D}}) where {B, D}
    return (1-d.dist.proba_zero) * d.counts
end

function get_proportion_observed(d::Dist)
    return d.counts
end

# function fit(zd::ZeroInflated, x::SampleZI)
#     if x.iszero
#         return ZeroInflated(zero(zd.dist), 1.0)
#     else
#         return ZeroInflated(fit(zd.dist, x.value), 0.0)
#     end
# end

# function fit(zd::ZeroInflated{D, F}, x) where {D, F}
#     if iszero(x)
#         return ZeroInflated(zero(zd.dist), 1.0)
#     else
#         return ZeroInflated(fit(zd.dist, x), 0.0)
#     end
# end
