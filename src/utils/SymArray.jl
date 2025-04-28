module FastSymArray

    export SymArray

    mutable struct SymArray{F}
        d::Dict{Tuple{Int, Int}, F}
        k::Int
    end

    function SymArray(k, d::F) where {F}
        @assert k > 0
        return SymArray{F}(Dict{Tuple{Int, Int}, F}(minmax(i, j) => d for i in 1:k
        for j in i:k), k)
    end

    Base.@propagate_inbounds function Base.getindex(a::SymArray, i, j)
        return a.d[minmax(i, j)]
    end

    function Base.setindex!(a::SymArray, v, i, j)
        a.d[minmax(i, j)] = v
    end


    function Base.sum(a::SymArray)
        return sum(values(a.d))
    end
end
