struct CategoricalSuffStats{M, T} <: SuffStats
    h::SVector{M, T}
    n::Int
end

function CategoricalSuffStats(num_categories::Int)
    h = SVector{num_categories, Int}(zeros(Int, num_categories))
    return CategoricalSuffStats{num_categories, Int}(h, 0)
end

function add_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] += 1
    ss = @set ss.n += 1
    return ss
end

function add_sample(ss::CategoricalSuffStats, ::Nothing)
    @reset ss.n += 1
    return ss
end

function remove_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] -= 1
    ss = @set ss.n -= 1
    return ss
end

function remove_sample(ss::CategoricalSuffStats, ::Nothing)
    @reset ss.n -= 1
    return ss
end

function make_k_block(k, ::Val{:categorical}; num_categories, kwargs...)
    k_block = SymArray{CategoricalSuffStats{num_categories, Int}}(undef, k, k)
    fill!(k_block, CategoricalSuffStats(num_categories))
    return k_block
end

function score(ss::CategoricalSuffStats)
    return ss.n - sum(abs2, ss.h) / max(ss.n, 1)
end

function to_params(ss::CategoricalSuffStats)
    n = max(ss.n, 1)
    return ss.h ./ n
end
