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


function check_compatiblity!(node_labels,g::GroupSize,)
    counts = StatsBase.countmap(node_labels)
    if length(counts) != g.number_groups || size(node_labels, 1) != sum(g)
        throw(ArgumentError("The vector of node labels is not compatible with the \
        group size: $(length(counts)) != $(g.number_groups) or $(size(node_labels, 1)) \
        != $(sum(g))"))
    end
    unbalanced = any(((k,v),) -> v != g[k], counts)
    if unbalanced
        @info "The group size is unbalanced, trying to fix it"
        g, node_labels = try_fixing_group_size!(node_labels, g)
        if any(((k,v),) -> v != g[k], StatsBase.countmap(node_labels))
            throw(ArgumentError("Could not fix the group size"))
        else
            @info "Fixed the group size by moving nodes between groups"
        end
    end
end


function try_fixing_group_size!(node_labels,g::GroupSize)
    counts = StatsBase.countmap(node_labels)
    groups_too_small = filter(((k,v),)  -> v < g[k], counts)
    groups_too_large = filter(((k,v),)  -> v > g[k], counts)
    amount_too_small = sum(g[k] - v for (k, v) in groups_too_small)
    amount_too_large = sum(v - g[k] for (k, v) in groups_too_large)
    if amount_too_small == amount_too_large
        nodes_to_move = []
        for (l,v) in groups_too_large
            number_nodes_to_move = v - g[l]
            nodes_to_move = vcat(nodes_to_move, findall(x-> x == l,node_labels)[1:number_nodes_to_move])
        end
        for (k, v) in groups_too_small
            number_nodes_to_move = g[k] - v
            for i in 1:number_nodes_to_move
                index = popfirst!(nodes_to_move)
                node_labels[index] = k
            end
        end
    end
    return g, node_labels
end
