abstract type AdditionalData end

struct BasicAdditionalInformation <: AdditionalData
    test::Int
end


struct Assignment{T} <: AbstractVector{Vector{Int}}
    group_size::GroupSize{T}
    node_labels::Vector{Int}
    additional_data::BasicAdditionalInformation
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

function get_edge_indices(a::Assignment,i::Int)
    nodes_i = get_vertex_in_group(a, i)
    return [(x, y) for x in nodes_i for y in nodes_i if x < y]
end


Base.size(a::Assignment) = (number_groups(a),number_groups(a))
Base.@propagate_inbounds function Base.getindex(a::Assignment, i::Int)
    @boundscheck checkbounds(a, i)
    return get_vertex_in_group(a, i)
end





function get_obs(g::SimpleGraph{T}, x::Tuple) where {T}
    return get_obs(g, x[1], x[2])
end

function get_obs(g::SimpleGraph{T}, i::Int, j::Int) where {T}
    return convert(Bool, has_edge(g, i, j))
end


function loglikelihood(a::Assignment, dists::SBM, g)
    loglikelihood = 0.0
    for i in 1:number_nodes(a)
        label_a = a.node_labels[i]
        for j in i+1:number_nodes(a)
            label_b = a.node_labels[j]
            loglikelihood += logdensityof(dists[label_a,label_b], get_obs(g,i,j))
        end
    end
    return loglikelihood
end


function fit(a::Assignment,g, distribution)
    sizes = [a.group_size[i] for i in 1:number_groups(a)]/number_nodes(a)
    dists = initialize_sbm(sizes, distribution)
    for group1 in 1:number_groups(a)
        for group2 in group1:number_groups(a)
            edges = get_edge_indices(a,group1,group2)
            dists[group1,group2] = fit(distribution, g, edges)
        end
    end
    return dists
end

function fit(distribution, g, edges)
    return Distributions.fit(typeof(distribution), get_obs.(Ref(g), edges))
end
