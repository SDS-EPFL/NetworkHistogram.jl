struct EdgeList{E}
    data::Vector{Vector{Tuple{Int,E}}}
end

function neighbors(A::EdgeList{E}, i::Int) where {E}
    return first.(A.data[i]), last.(A.data[i])
end

function iterate_neighbors(A::EdgeList{E}, i::Int) where {E}
    return zip(first.(A.data[i]), last.(A.data[i]))
end

function edge_type(edgelist::EdgeList{E}) where {E}
    return E
end

function nodes(edgelist::EdgeList{E}) where {E}
    return length(edgelist.data)
end


function EdgeList(A::AbstractMatrix{<:Union{Missing,E}}) where {E}
    n = size(A, 1)
    data = Vector{Vector{Tuple{Int,E}}}(undef, n)
    for j in 1:n
        data[j] = Vector{Tuple{Int,E}}(undef, 0)
        for i in 1:n
            if !ismissing(A[i,j])
                push!(data[j], (i, A[i, j]))
            end
        end
    end
    return EdgeList(data)
end


function Base.convert(::Type{EdgeList{E}}, A::AbstractMatrix{E}) where {E}
    return EdgeList(A)
end


function fit(d::Dist, A::EdgeList{E}) where {E}
    new_data = Vector{Vector{Tuple{Int, typeof(d)}}}(undef, length(A.data))
    for j in 1:length(A.data)
        new_data[j] = Vector{Tuple{Int, typeof(d)}}(undef, length(A.data[j]))
        for (k,(i, e)) in enumerate(A.data[j])
            new_data[j][k] = (i, fit(d, e))
        end
    end
    return EdgeList(new_data)
end
