"""
Array-like storage for the number of nodes in each group. Try to split the number of nodes
into equal groups, but if it is not possible, the last group may have more nodes.
"""
struct GroupSize{T} <: AbstractVector{Int}
    group_number::T
    number_groups::Int

    function GroupSize(number_nodes, h::Real)
        @assert 0 < h < 1
        standard_group = floor(Int, number_nodes * h)
        GroupSize(number_nodes, standard_group)
    end

    function GroupSize(number_nodes, standard_group::Integer)
        @assert 1 < standard_group <= number_nodes
        number_groups = number_nodes ÷ standard_group # number of standard groups!
        if number_groups * standard_group == number_nodes
            new{Int}(standard_group, number_groups)
        else
            remainder_group = standard_group +
                              mod(number_nodes, standard_group)
            new{Tuple{Int, Int}}(
                (standard_group, remainder_group), number_groups)
        end
    end
end

Base.size(g::GroupSize) = (g.number_groups,)
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Int}, i::Int)
    @boundscheck checkbounds(g, i)
    return g.group_number
end

Base.@propagate_inbounds function Base.getindex(
        g::GroupSize{Tuple{Int, Int}}, i::Int)
    @boundscheck checkbounds(g, i)
    return i < length(g) ? g.group_number[1] : g.group_number[2]
end

mutable struct Assignment{E, D, F}
    node_labels::AbstractVector{Int}
    const edges::EdgeList{E}
    const dists::EdgeList{D}
    θ::SymArray{D}
    log_likelihood::SymArray{F}
end

number_nodes(a::Assignment) = length(a.node_labels)
number_groups(a::Assignment) = size(a.θ, 1)

function loglikelihood(a::Assignment)
    return sum(a.log_likelihood)
end

function group(a::Assignment, node::Int)
    return a.node_labels[node]
end

function get_edges_in_groups(a::Assignment, g1::Int, g2::Int)
    return get_edges_in_groups(a.node_labels, a.edges, g1, g2)
end

function get_edges_in_groups(node_labels, edges_all, g1, g2)

    edges = Vector{edge_type(edges_all)}()
    nodes_g1 = findall(x -> x == g1, node_labels)
    nodes_g2 = findall(x -> x == g2, node_labels)

    for u in nodes_g1
        for (v, e) in iterate_neighbors(edges_all, u)
            if  v in nodes_g2 && ((g1 == g2 && u < v) || g1 != g2)
                push!(edges, e)
            end
        end
    end
    return edges
end

function Assignment(node_labels, edge_list::EdgeList{E}, dist::Dist{D}) where {E, D}
    dists = fit(dist, edge_list)
    number_groups = length(unique(node_labels))
    θ = SymArray(number_groups, dist)
    log_likelihood = SymArray(number_groups, 0.0)
    for u in 1:nodes(dists)
        g1 = node_labels[u]
        for (v, d) in iterate_neighbors(dists, u)
            g2 = node_labels[v]
            if u < v
                θ[g1, g2] = add_to(θ[g1, g2], d)
            end
        end
    end
    for k in 1:number_groups
        for l in k:number_groups
            log_likelihood[k,
                l] = loglikelihood(
                θ[k, l], get_edges_in_groups(node_labels, edge_list, k, l))
        end
    end
    return Assignment(node_labels, edge_list, dists, θ, log_likelihood)
end
