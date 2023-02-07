"""
Array-like storage for the number of nodes in each group.
"""
struct GroupSize{T} <: AbstractVector{Int}
    group_number::T
    number_groups::Int

    function GroupSize(number_nodes, h::Real)
        @assert 0 < h < 1
        standard_group = floor(Int, number_nodes * h)
        GroupSize(number_nodes, standard_group)
    end

    function GroupSize(number_nodes, standard_group::Int)
        @assert 1 < standard_group < number_nodes
        number_groups = number_nodes ÷ standard_group # number of standard groups!
        if number_groups * standard_group == number_nodes
            new{Int}(standard_group, number_groups)
        else
            remainder_group = number_nodes - number_groups * standard_group
            if remainder_group == 1
                @warn "h has to be changes as only one node in remainder group"
                standard_group -= 1
                remainder_group = number_groups + 1 # because equal to 1+number_groups because we take 1 from each standard group, and there are number_groups of them
                if standard_group == 1
                    error("Standard group size now 1, please choose a new h value.")
                end
            end
            new{Tuple{Int, Int}}((standard_group, remainder_group), number_groups + 1)
        end
    end
end

Base.size(g::GroupSize) = (g.number_groups,)
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Int}, i::Int)
    @boundscheck checkbounds(g, i)
    return g.group_number
end

Base.@propagate_inbounds function Base.getindex(g::GroupSize{Tuple{Int, Int}}, i::Int)
    @boundscheck checkbounds(g, i)
    return i < length(g) ? g.group_number[1] : g.group_number[2]
end
