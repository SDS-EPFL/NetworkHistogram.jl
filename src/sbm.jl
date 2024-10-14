struct SBM{T,K} <: AbstractMatrix{T}
    sizes :: Vector{Float64}
    probs:: SymmetricTensor{T, K, K}

end



function initialize_sbm(sizes::Vector{Float64}, dist, k = length(sizes))
    probs = Vector{typeof(dist)}(undef, binomial(k + 1, 2))
    fill!(probs, dist)
    return SBM(sizes, SymmetricTensor(probs, Val(k), Val(k)))
end

function initialize_sbm(k::Int, dist)
    return initialize_sbm(ones(k) / k, dist)
end

number_blocks(::SBM{T,K}) where {T,K} = K

Base.size(s::SBM)= size(s.probs)
Base.ndims(s::SBM) = 2
Base.eltype(s::SBM{T,K}) where {T,K} = T
Base.setindex!(s::SBM, v, i, j) = setindex!(s.probs, v, i, j)
Base.@propagate_inbounds function Base.getindex(s::SBM, i, j)
    @boundscheck checkbounds(s.probs, i, j)
    return Base.getindex(s.probs, i, j)
end
