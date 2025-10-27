
abstract type GenericSuffStatsType <: SuffStats end

struct GenericSuffStats{T} <: GenericSuffStatsType
    samples::Vector{T}
end

function GenericSuffStats(::AbstractArray{T}) where {T}
    return GenericSuffStats{T}(Vector{T}())
end

function get_samples(ss::GenericSuffStats)
    return ss.samples
end

function add_sample(ss::GenericSuffStats, sample)
    append!(ss.samples, sample)
    return ss
end

function remove_sample(ss::GenericSuffStats, sample)
    index = findfirst(==(sample), ss.samples)
    if index !== nothing
        deleteat!(ss.samples, index)
    end
    return ss
end

function make_k_block(k, generic; data::AbstractArray, kwargs...)
    @warn "Using GenericSuffStats may lead to high memory usage for large datasets.
         Consider using more specialized sufficient statistics types when possible."
    k_block = SymArray{GenericSuffStats{eltype(data)}}(undef, k, k)
    for j in 1:k, i in 1:k
        k_block[i, j] = GenericSuffStats(data)
    end
    return k_block
end

# use indices rather than pushing and deleting samples for better performance ?
# struct GenericSuffStatsIndex{T} <: GenericSuffStatsType
#     indices::Vector{Tuple{Int, Int}}
#     data::T
# end

# function get_samples(ss::GenericSuffStatsIndex)
#     return [ss.data[i, j] for (i, j) in ss.indices]
# end

# function GenericSuffStatsIndex{T}(data::T) where {T}
#     return GenericSuffStatsIndex{T}(Vector{Tuple{Int, Int}}(), data)
# end

# function add_sample(ss::GenericSuffStatsIndex, sample, i, j)
#     push!(ss.indices, (i, j))
#     return ss
# end

function score(ss::GenericSuffStatsType; dist::D, kwargs...) where {D}
    if dist === nothing
        @error("No distribution provided for scoring GenericSuffStats")
    end
    samples = get_samples(ss)
    d = fit(D, samples)
    return -sum(logpdf.(d, samples))
end
