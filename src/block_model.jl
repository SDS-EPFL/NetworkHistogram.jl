struct BlockModel{D, V}
    _dists::SymArray{D}
    sizes::V
    cum_sizes::V
end

function BlockModel(k::Int, d::D) where {D}
    sizes = fill(1 / k, k)
    cumulative_sizes = cumsum(sizes)
    _dists = SymArray(k, d)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

function BlockModel(a::Assignment)
    k = length(unique(a.node_labels))
    sizes = proportions(a)
    cumulative_sizes = cumsum(sizes)
    _dists = unwrap.(a.θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

function BlockModel(nodes_labels, θ)
    k = length(unique(nodes_labels))
    sizes = counts(nodes_labels) / length(nodes_labels)
    cumulative_sizes = cumsum(sizes)
    _dists = unwrap.(θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

function BlockModel(θ::AbstractMatrix{D}) where {D}
    k = size(θ, 1)
    sizes = fill(1 / k, k)
    cumulative_sizes = cumsum(sizes)
    _dists = convert(SymArray{D}, θ)
    return BlockModel(_dists, sizes, cumulative_sizes)
end

function map_ξ_to_block(bm::BlockModel, ξ::T) where {T <: Real}
    return findfirst(x -> x >= ξ, bm.cum_sizes)
end

function sample(bm::BlockModel, latents::Int, args...)
    latents = map(x -> map_ξ_to_block(bm, x), rand(latents))
    return latents, sample(bm, latents, args...)
end

function sample(bm::BlockModel, latents::Vector{T}, args...) where {T}
    A = Array{eltype(bm[1, 1]), 2}(undef, length(latents), length(latents))
    for j in 1:length(latents)
        for i in 1:(j - 1)
            A[i, j] = A[j, i]
        end
        for i in (j + 1):length(latents)
            # println("latents[i]: ", latents[i], " latents[j]: ", latents[j])
            # println("bm[latents[i], latents[j]]: ", bm[latents[i], latents[j]])
            A[i, j] = sample(bm[latents[i], latents[j]], args...)
            A[j, i] = A[i, j]
        end
    end
    # fill the diagonal with zeros, avoid undefined references
    for i in 1:length(latents)
        A[i, i] = zero(A[1, 2])
    end
    return A
end

# this is probably awfull

function Base.getindex(s::BlockModel, i::Int, j::Int)
    return s._dists[i, j]
end

function Base.setindex!(s::BlockModel, v, i::Int, j::Int)
    s._dists[i, j] = v
end

function Base.size(s::BlockModel)
    return (s._dists.k, s._dists.k)
end

function Base.getindex(s::BlockModel, i::Real, j::Real)
    k = findfirst(x -> x ≥ i, s.cum_sizes)
    l = findfirst(x -> x ≥ j, s.cum_sizes)
    return s._dists[k, l]
end

function Base.setindex!(s::BlockModel, v, i::Real, j::Real)
    k = findfirst(x -> x ≥ i, s.cum_sizes)
    l = findfirst(x -> x ≥ j, s.cum_sizes)
    s._dists[k, l] = v
end

# helpers for generating ordered latents
function ordered_latents(bm::BlockModel, n::Int)
    return sort(map(x -> map_ξ_to_block(bm, x), rand(n)))
end

function get_probability_matrix(
        bm::BlockModel{D}, latents::AbstractVector, default_dist = nothing) where {D}
    # hack for dirac at 0 dist (no self-loop)
    if isnothing(default_dist)
        try
            default_dist = zero(bm[1, 1])
        catch e
            if !is(e, MethodError)
                rethrow(e)
            end
            error("Please provide a default distribution for the diagonal as it could not be inferred")
        end
    end
    n = length(latents)
    A = Array{D, 2}(undef, n, n)
    for j in 1:n
        for i in 1:n
            if i == j
                A[i, i] = default_dist
            else
                A[i, j] = bm[latents[i], latents[j]]
            end
        end
    end
    return A
end

function get_probability_matrix(a::Assignment, default_dist = nothing)
    return get_probability_matrix(BlockModel(a.θ), a.node_labels, default_dist)
end
