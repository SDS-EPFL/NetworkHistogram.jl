include("group_numbering.jl")

"""
    struct Assignment{T, B} <: AbstractVector{Vector{Int}}

A structure representing an assignment of nodes to groups.

# Fields
- `group_size::GroupSize{T}`: The size of each group.
- `node_labels::Vector{Int}`: A vector of node labels.
- `additional_data::B`: Additional data associated with the assignment.

# Constructor
    Assignment(group_size::GroupSize{T}, node_labels, additional_data::B) where {T, B}

Creates a new `Assignment` instance.

# Arguments
- `group_size::GroupSize{T}`: The size of each group.
- `node_labels::Vector{Int}`: A vector of node labels. The length of this vector must be equal to the sum of `group_size`.
- `additional_data::B`: Additional data associated with the assignment.

# Throws
- `ArgumentError`: If the length of `node_labels` is not equal to the sum of `group_size`.
"""
struct Assignment{T, B} <: AbstractVector{Vector{Int}}
    group_size::GroupSize{T}
    node_labels::Vector{Int}
    additional_data::B

    function Assignment(group_size::GroupSize{T}, node_labels,
            additional_data::B) where {T, B}
        if length(node_labels) != sum(group_size)
            throw(ArgumentError("The length of `node_labels` must be equal to the sum of \
                                `group_size`"))
        end
        return new{T, B}(group_size, node_labels, additional_data)
    end
end

function Assignment(group_size::GroupSize, node_labels)
    if length(node_labels) != sum(group_size)
        throw(ArgumentError("The length of `node_labels` $(length(node_labels)) must be \
        equal to the sum of `group_size` $(sum(group_size))"))
    end
    c = StatsBase.countmap(node_labels)
    if length(c) != length(group_size)
        throw(ArgumentError("The number of unique elements in `node_labels` $(length(c)) \
        must be equal to the length of `group_size` $(length(group_size))"))
    end
    for (k, v) in c
        if v != group_size[k]
            throw(ArgumentError("The number of elements in `node_labels` $(v) for group \
            $(k) must be equal to the size of the group $(group_size[k])"))
        end
    end
    return Assignment(group_size, node_labels, nothing)
end

function number_groups(assignment::Assignment)
    return length(assignment.group_size)
end

function number_nodes(assignment::Assignment)
    return length(assignment.node_labels)
end

function get_vertex_in_group(assignment::Assignment, group)
    return findall(assignment.node_labels .== group)
end

function get_group_of_vertex(assignment::Assignment, vertex)
    return assignment.node_labels[vertex]
end

function get_edge_indices(a::Assignment, i, j)
    if i == j
        return get_edge_indices(a, i)
    else
        return [(x, y) for x in get_vertex_in_group(a, i)
                for y in get_vertex_in_group(a, j)]
    end
end

function get_edge_indices(a::Assignment, i)
    nodes_i = get_vertex_in_group(a, i)
    return [(x, y) for x in nodes_i for y in nodes_i if x < y]
end

Base.size(a::Assignment) = (number_groups(a),)
Base.@propagate_inbounds function Base.getindex(a::Assignment, i)
    @boundscheck checkbounds(a, i)
    return get_vertex_in_group(a, i)
end

function get_ordered_adjacency_matrix(a::Assignment)
    perm = sortperm(a.node_labels)
    return a.additional_data.A[perm, perm]
end
