struct BlockModel{D, K, T}
    _dists::SymArray{D}
    sizes::SVector{K, T}
    cum_sizes::SVector{K,T}
end

function BlockModel(k::Int, d::D) where {D}
    sizes = @SVector fill(1/k, k)
    cumulative_sizes = SVector{k}(cumsum(sizes))
    _dists = SymArray(k, d)
    return BlockModel{D, k, Float64}(_dists, sizes, cumulative_sizes)
end


function BlockModel(a::Assignment)
    k = length(unique(a.node_labels))
    sizes = SVector{k}(proportions(a))
    cumulative_sizes = SVector{k}(cumsum(sizes))
    _dists = unwrap.(a.θ)
    return BlockModel{eltype(_dists), k, eltype(cumulative_sizes)}(_dists, sizes, cumulative_sizes)
end


function sample(bm::BlockModel, latents::Vector{T}) where {T}
    A = Array{eltype(bm[1,1]), 2}(undef, length(latents), length(latents)) .* zero(eltype(bm[1,1]))
    for j in 1:length(latents)
        for i in 1:j-1
            A[i, j] = A[j, i]
        end
        for i in j+1:length(latents)
            # println("i: ", i, " j: ", j)
            # println("latents[i]: ", latents[i], " latents[j]: ", latents[j])
            # println("bm[latents[i], latents[j]]: ", bm[latents[i], latents[j]])
                A[i, j] = sample(bm[latents[i], latents[j]])
                A[j, i] = A[i, j]
        end
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
