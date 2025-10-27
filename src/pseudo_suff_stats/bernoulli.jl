
struct BernoulliSuffStats{T} <: SuffStats
    h::T
    n::T
end

function BernoulliSuffStats()
    return BernoulliSuffStats{Int}(0, 0)
end

@inline function add_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h += 1)
    @reset ss.n += 1
    return ss
end

@inline function add_sample(ss::BernoulliSuffStats, ::Nothing)
    @reset ss.n += 1
    return ss
end

@inline function remove_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h -= 1)
    @reset ss.n -= 1
    return ss
end

@inline function remove_sample(ss::BernoulliSuffStats, ::Nothing)
    @reset ss.n -= 1
    return ss
end

function make_k_block(k, ::Val{:binary}; kwargs...)
    k_block = SymArray{BernoulliSuffStats{Int}}(undef, k, k)
    fill!(k_block, BernoulliSuffStats())
    return k_block
end

@inline function score(ss::BernoulliSuffStats)
    n = max(ss.n, 1)
    p = ss.h / n
    return n * (xlogx(1 - p) + xlogx(p))
end
