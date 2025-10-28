struct CategoricalSuffStats{M, T} <: SuffStats
    h::SVector{M, T}
end

function CategoricalSuffStats(num_categories::Int)
    h = SVector{num_categories, Int}(zeros(Int, num_categories))
    return CategoricalSuffStats{num_categories, Int}(h)
end

function add_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] += 1
    return ss
end

function remove_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] -= 1
    return ss
end

function make_k_block(k, ::Val{:categorical}; num_categories, kwargs...)
    k_block = SymArray{CategoricalSuffStats{num_categories, Int}}(undef, k, k)
    fill!(k_block, CategoricalSuffStats(num_categories))
    return k_block
end

function score(ss::CategoricalSuffStats)
    n = sum(ss.h)
    return n - sum(abs2, ss.h) / max(n, 1)
end

function to_params(ss::CategoricalSuffStats)
    return custom_normalize(ss.h)
end

function custom_normalize(ps::SVector{M, T}) where {M, T}
    n = sum(ps)
    n == 0 && return zero(SVector{M, T})
    return ps / n
end
