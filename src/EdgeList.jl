"""
    EdgeList{E}

A memory-efficient adjacency list representation for sparse networks.

# Fields
- `data::Vector{Vector{E}}`: For each node, stores the edge values to its neighbors
- `name_list::Vector{Vector{Int}}`: For each node, stores the node indices of its neighbors

# Type Parameters
- `E`: The type of edge values (e.g., Int, Float64, or custom distribution types)

# Examples
```julia
# From an adjacency matrix
A = [0 1 0; 1 0 1; 0 1 0]
edges = EdgeList(A)

# Access neighbors of node 1
neighbor_indices, edge_values = neighbors(edges, 1)

# Iterate through neighbors
for (neighbor, edge) in iterate_neighbors(edges, 1)
    println("Edge to node ", neighbor, " with value ", edge)
end
```

See also: [`neighbors`](@ref), [`iterate_neighbors`](@ref), [`get_edge`](@ref)
"""
struct EdgeList{E}
    data::Vector{Vector{E}}
    name_list::Vector{Vector{Int}}
end

"""
    neighbors(A::EdgeList, i::Int)

Get the neighbor indices and edge values for node `i`.

Returns a tuple `(neighbor_indices, edge_values)` where each vector has the same length.

# Example
```julia
edges = EdgeList(A)
neighbor_nodes, edge_vals = neighbors(edges, 1)
```
"""
@inline function neighbors(A::EdgeList{E}, i::Int) where {E}
    @boundscheck checkbounds(A.data, i)
    @boundscheck checkbounds(A.name_list, i)
    return A.name_list[i], A.data[i]
end

"""
    iterate_neighbors(A::EdgeList, i::Int)

Returns an iterator over (neighbor_index, edge_value) pairs for node `i`.

# Example
```julia
for (j, edge) in iterate_neighbors(edges, i)
    # Process edge from i to j
end
```
"""
@inline iterate_neighbors(A::EdgeList, i::Int) = zip(neighbors(A, i)...)

"""
    edge_type(A::EdgeList{E})

Get the element type `E` of edges stored in the EdgeList.
"""
@inline edge_type(A::EdgeList{E}) where {E} = E

"""
    nodes(edgelist::EdgeList)
    number_nodes(edgelist::EdgeList)

Return the number of nodes in the network.
"""
@inline nodes(edgelist::EdgeList) = length(edgelist.data)
@inline number_nodes(edgelist::EdgeList) = nodes(edgelist)

"""
    EdgeList(A::AbstractMatrix{<:Union{Missing, E}}) where {E}

Construct an EdgeList from an adjacency matrix. Missing values are treated as absent edges,
and diagonal entries are excluded (no self-loops).

# Arguments
- `A::AbstractMatrix`: Adjacency matrix where `missing` indicates absent edges

# Example
```julia
A = [0 1 missing; 1 0 2; missing 2 0]
edges = EdgeList(A)
```
"""
function EdgeList(A::AbstractMatrix{<:Union{Missing, E}}) where {E}
    _from_adj_to_edge_list(A)
end
EdgeList(adj_list::EdgeList) = adj_list

"""
    get_edge(A::EdgeList{E}, i::Int, j::Int) where {E}

Get the edge value between nodes `i` and `j`. Returns `zero(E)` if no edge exists or if `i == j`.

# Arguments
- `A::EdgeList{E}`: The edge list
- `i::Int`: Source node index
- `j::Int`: Target node index

# Returns
- Edge value of type `E`, or `zero(E)` if no edge exists
"""
function get_edge(A::EdgeList{E}, i::Int, j::Int) where {E}
    if i == j
        return zero(E)
    end
    if j ∉ A.name_list[i] && i ∉ A.name_list[j]
        return zero(E)
    end
    for (k, e) in iterate_neighbors(A, i)
        if k == j
            return e
        end
    end
    return zero(E)  # If edge not found in the iteration
end

# Internal function to convert adjacency matrix to EdgeList format
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
            # Exclude diagonal and missing edges
            if !ismissing(A[i, j]) && i != j
                push!(name_list[j], i)
                push!(data[j], function_to_apply(A[i, j]))
            end
        end
    end
    return EdgeList(data, name_list)
end

# Internal functions for preprocessing edge data
function _fast_compressed_obs(d::Dist, A::AbstractMatrix, zeroinflated)
    _from_adj_to_edge_list(A, x -> _fast_compressed_obs(d, x, zeroinflated))
end
function _fast_compressed_obs(d::Dist, A::EdgeList{E}, zeroinflated) where {E}
    _make_shift_broadcast(A.data, x -> _fast_compressed_obs(d, x, zeroinflated))
end

# Internal function to apply a transformation to EdgeList data
function _make_shift_broadcast(A::EdgeList, f)
    n = length(A.data)
    test = f(A.data[1][1])
    data = Vector{Vector{typeof(test)}}(undef, n)
    for j in 1:n
        data[j] = f.(A.data[j])
    end
    return EdgeList(data, A.name_list)
end

"""
    fit(d::Dist, A::EdgeList{E}) where {E}

Fit the distribution `d` to each edge in the EdgeList `A`, returning a new EdgeList
where each edge is replaced by its fitted distribution.

# Arguments
- `d::Dist`: The distribution type to fit
- `A::EdgeList{E}`: EdgeList containing edge observations

# Returns
- `EdgeList{typeof(d)}`: New EdgeList with fitted distributions
"""
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
