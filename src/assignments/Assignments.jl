struct Assignment{T, B} <: AbstractVector{Vector{Int}}
    group_size::GroupSize{T}
    node_labels::Vector{Int}
    additional_data::B

    function Assignment(group_size::GroupSize{T}, node_labels::Vector{Int},
            additional_data::B) where {T, B}
        if length(node_labels) != sum(group_size)
            throw(ArgumentError("The length of `node_labels` must be equal to the sum of `group_size`"))
        end
        return new{T, B}(group_size, node_labels, additional_data)
    end
end

function Assignment(group_size::GroupSize{T}, node_labels::Vector{Int}) where {T}
    if length(node_labels) != sum(group_size)
        throw(ArgumentError("The length of `node_labels` must be equal to the sum of `group_size`"))
    end
    return Assignment(group_size, node_labels, nothing)
end

function number_groups(assignment::Assignment)
    return length(assignment.group_size)
end

function number_nodes(assignment::Assignment)
    return length(assignment.node_labels)
end

function get_vertex_in_group(assignment::Assignment, group::Int)
    return findall(assignment.node_labels .== group)
end

function get_edge_indices(a::Assignment, i::Int, j::Int)
    return [(x, y) for x in get_vertex_in_group(a, i)
            for y in get_vertex_in_group(a, j)]
end

function get_edge_indices(a::Assignment, i::Int)
    nodes_i = get_vertex_in_group(a, i)
    return [(x, y) for x in nodes_i for y in nodes_i if x < y]
end

Base.size(a::Assignment) = (number_groups(a),)
Base.@propagate_inbounds function Base.getindex(a::Assignment, i::Int)
    @boundscheck checkbounds(a, i)
    return get_vertex_in_group(a, i)
end

function loglikelihood(a::Assignment, dists::SBM, g)
    loglikelihood = 0.0
    for i in 1:number_nodes(a)
        label_a = a.node_labels[i]
        for j in (i + 1):number_nodes(a)
            label_b = a.node_labels[j]
            loglikelihood += logdensityof(dists[label_a, label_b], get_obs(g, i, j))
        end
    end
    return loglikelihood
end

function fit(a::Assignment, g, distribution)
    dists = initialize_sbm(a.group_size, distribution)
    for group1 in 1:number_groups(a)
        for group2 in group1:number_groups(a)
            edges = get_edge_indices(a, group1, group2)
            dists[group1, group2] = fit(distribution, g, edges)
        end
    end
    return dists
end

function fit(distribution, g, edges)
    return Distributions.fit(typeof(distribution), get_obs.(Ref(g), edges))
end
