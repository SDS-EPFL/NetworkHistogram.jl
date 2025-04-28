struct BlockModel{D,K,T}
    _dists::SymArray{D}
    sizes::SVector{K,T}
end


Base.@propagate_inbounds function Base.getindex(s::BlockModel, i, j)
    return s._dists[minmax(i, j)]
end

function Base.setindex!(s::BlockModel, v, i, j)
    s._dists[minmax(i, j)] = v
end
