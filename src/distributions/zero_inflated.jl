
struct ZeroInflated{D, F}
    dist::D
    proba_zero::F
end

abstract type SampleZI end

struct SampleNZ{F} <: SampleZI
    value::F
end

struct SampleZ <: SampleZI
end

function ZeroInflated(dist)
    return ZeroInflated(dist, 0.0)
end

function logpdf(zi::ZeroInflated{D, F}, x::SampleZ) where {D, F}
    return log(zi.proba_zero)
end

function logpdf(zi::ZeroInflated{D, F}, x::SampleNZ{T}) where {D, F, T}
    return log(1 - zi.proba_zero) + logpdf(zi.dist, x.value)
end

function agg_params(
        zi1::ZeroInflated{D, F}, zi2::ZeroInflated{D, F}, w1, w2) where {D, F}
    new_dist = agg_params(zi1.dist, zi2.dist, w1, w2)
    new_proba_zero = w1 * zi1.proba_zero + w2 * zi2.proba_zero
    return ZeroInflated(new_dist, new_proba_zero)
end

zero(zi::ZeroInflated) = ZeroInflated(zero(zi.dist), 0.0)

eltype(zi::ZeroInflated{D, F}) where {D, F} = Union{SampleZ, SampleNZ{eltype(D)}}


function fit(zi::ZeroInflated{D, F}, x::SampleZ) where {D, F}
    return ZeroInflated(zero(zi.dist), 1.0)
end

function fit(zi::ZeroInflated{D, F}, x::SampleNZ{T}) where {D, F, T}
    return ZeroInflated(fit(zi.dist, x.value), 0.0)
end

function _fast_compressed_obs(zi::ZeroInflated, x)
    if iszero(x)
        return SampleZ()
    else
        return SampleNZ(_fast_compressed_obs(zi.dist, x))
    end
end
