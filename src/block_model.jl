struct BlockModel{D, K, T}
    _dists::SymArray{D}
    sizes::SVector{K, T}
    cum_sizes::Vector{T}
end

function BlockModel(k::Int, d::D) where {D}
    sizes = @SVector fill(1/k, k)
    cumulative_sizes = cumsum(sizes)
    _dists = SymArray(k, d)
    return BlockModel{D, k, Float64}(_dists, sizes, cumulative_sizes)
end


function sample(bm::BlockModel, latents::Vector{T}) where {T}
    #fuck need the element type of the distribution...
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
