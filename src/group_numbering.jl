struct GroupSize{T} <: AbstractVector{Int}
    group_number::T
    num_groups::Int
    function GroupNumber(number_nodes, h::Float)
        @assert 0 < h < 1
        standard_group = floor(Int, number_nodes * h)
        GroupNumber(number_nodes, standard_group)
    end
    function GroupNumber(number_nodes, standard_group::Int)
        @assert 0 < standard_group < number_nodes
        num_groups = number_nodes ÷ standard_group
        if number_groups * standard_group == number_nodes
            new{Int}(standard_group, num_groups)
        else
            extra_group = number_nodes - number_groups * standard_group
            new{Tuple{Int,Int}}((standard_group, extra_group), num_groups+1)
        end
    end
end
Base.size(g::GroupSize) = (g.num_groups,)
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Int}, i::Int)
    @boundscheck checkbounds(g, i)
    return g.group_number
end
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Tuple{Int,Int}}, i::Int)
    @boundscheck checkbounds(g, i)
    return i < length(g) ? g.group_number[1] : g.group_number[2]
end