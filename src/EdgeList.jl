struct EdgeList{E}
    data::Vector{Vector{E}}
    name_list::Vector{Vector{Int}}
end

function neighbors(A::EdgeList{E}, i::Int) where {E}
    return A.name_list[i], A.data[i]
end

iterate_neighbors(A::EdgeList, i::Int) = zip(neighbors(A, i)...)
edge_type(A::EdgeList{E}) where {E} = E
nodes(edgelist::EdgeList) = length(edgelist.data)
number_nodes(edgelist::EdgeList) = nodes(edgelist)

function EdgeList(A::AbstractMatrix{<:Union{Missing, E}}) where {E}
    _from_adj_to_edge_list(A)
end
EdgeList(adj_list::EdgeList) = adj_list

function get_edge(A::EdgeList{E}, i::Int, j::Int) where {E}
    if i == j
        return zero(E)
    end
    # TODO: probably can remove this
    if j ∉ A.name_list[i] && i ∉ A.name_list[j]
        return zero(E)
    end
    for (k, e) in iterate_neighbors(A, i)
        if k == j
            return e
        end
    end
end

# function EdgeList(A::AbstractMatrix{<:Union{Missing,E}}) where {E}
#     n = size(A, 1)
#     data = Vector{Vector{E}}(undef, n)
#     name_list = Vector{Vector{Int}}(undef, n)
#     for j in 1:n
#         data[j] = Vector{E}(undef, 0)
#         name_list[j] = Vector{Int}(undef, 0)
#         for i in 1:n
#             if !ismissing(A[i,j]) # gonna be an issue with MC! have to define 0 chain and fast operations on them
#                 push!(name_list[j], i)
#                 push!(data[j], A[i, j])
#             end
#         end
#     end
#     return EdgeList(data, name_list)
# end

function _from_adj_to_edge_list(
        A::AbstractMatrix, function_to_apply = identity)
    n = size(A, 1)
    input = findfirst(x -> !ismissing(x), A)
    test = function_to_apply(A[input])
    data = Vector{Vector{typeof(test)}}(undef, n)
    name_list = Vector{Vector{Int}}(undef, n)
    for j in 1:n
        data[j] = Vector{typeof(test)}(undef, 0)
        name_list[j] = Vector{Int}(undef, 0)
        for i in 1:n
            if !ismissing(A[i, j])
            end
            if !ismissing(A[i, j]) && i != j # gonna be an issue with MC! have to define 0 chain and fast operations on them
                push!(name_list[j], i)
                push!(data[j], function_to_apply(A[i, j]))
            end
        end
    end
    return EdgeList(data, name_list)
end

function _fast_compressed_obs(d::Dist, A::AbstractMatrix)
    _from_adj_to_edge_list(A, x -> _fast_compressed_obs(d, x))
end
function _fast_compressed_obs(d::Dist, A::EdgeList{E}) where {E}
    _make_shift_broadcast(A.data, x -> _fast_compressed_obs(d, x))
end

function _make_shift_broadcast(A::EdgeList, f)
    # may work ? -> data = f.(A.data)
    n = length(A.data)
    test = f(A.data[1][1])
    data = Vector{Vector{typeof(test)}}(undef, n)
    for j in 1:n
        data[j] = f.(A.data[j])
    end
    return EdgeList(data, A.name_list)
end

#convert(::Type{EdgeList}, A::AbstractMatrix) =  EdgeList(A)

function fit(d::Dist, A::EdgeList{E}) where {E}
    new_data = Vector{Vector{typeof(d)}}(undef, length(A.data))
    for j in 1:length(A.data)
        new_data[j] = Vector{typeof(d)}(undef, length(A.data[j]))
        for (k, e) in enumerate(A.data[j])
            new_data[j][k] = fit(d, e)
        end
    end
    return EdgeList(new_data, A.name_list)
end
