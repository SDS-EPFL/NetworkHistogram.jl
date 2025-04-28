struct EdgeList{E}
    data::Vector{Vector{Tuple{Int,E}}}
end

function neighbors(A::EdgeList{E}, i::Int) where {E}
    return first.(A.data[i]), last.(A.data[i])
end

function Base.eltype(edgelist::EdgeList{E}) where {E}
    return E
end

function nodes(edgelist::EdgeList{E}) where {E}
    return length(edgelist.data)
end


function EdgeList(A::AbstractMatrix{E}) where {E}
    n = size(A, 1)
    data = Vector{Vector{Tuple{Int,E}}}(undef, n)
    for j in 1:n
        data[j] = Vector{Tuple{Int,E}}(undef, 0)
        for i in 1:n
            if A[i, j] != 0
                push!(data[j], (i, A[i, j]))
            end
        end
    end
    return EdgeList(data)
end
