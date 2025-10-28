struct GenericSuffStats{T, D} <: SuffStats
    samples::Vector{T}
    dist::D
end

function GenericSuffStats(::AbstractArray{T}, dist::D) where {T, D}
    return GenericSuffStats{T, D}(Vector{T}(), dist)
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

function make_k_block(k, generic; data::AbstractArray, dist::D, kwargs...) where {D}
    @warn "Using GenericSuffStats may be very slow even for small graphs.
         Consider using more specialized sufficient statistics types when possible."
    k_block = SymArray{GenericSuffStats{eltype(data), D}}(undef, k, k)
    for j in 1:k, i in 1:k
        k_block[i, j] = GenericSuffStats(data, dist)
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

function loss(ss::GenericSuffStats{T, D}) where {T, D}
    samples = get_samples(ss)
    d = fit(D, samples)
    return -sum(logpdf.(d, samples))
end

function to_params(ss::GenericSuffStats)
    d = fit(typeof(ss.dist), get_samples(ss))
    return params(d)
end
