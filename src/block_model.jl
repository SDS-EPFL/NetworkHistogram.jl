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
