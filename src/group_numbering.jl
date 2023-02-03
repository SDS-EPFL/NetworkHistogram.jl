struct GroupSize{T} <: AbstractVector{Int}
    group_number::T
    number_groups::Int
    function GroupSize(number_nodes, h::Real)
        @assert 0 < h < 1
        standard_group = floor(Int, number_nodes * h)
        GroupSize(number_nodes, standard_group)
    end
    function GroupSize(number_nodes, standard_group::Int)
        @assert 0 < standard_group < number_nodes
        number_groups = number_nodes ÷ standard_group
        if number_groups * standard_group == number_nodes
            new{Int}(standard_group, number_groups)
        else
            extra_group = number_nodes - number_groups * standard_group
            new{Tuple{Int,Int}}((standard_group, extra_group), number_groups+1)
        end
    end
end
Base.size(g::GroupSize) = (g.number_groups,)
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Int}, i::Int)
    @boundscheck checkbounds(g, i)
    return g.group_number
end
Base.@propagate_inbounds function Base.getindex(g::GroupSize{Tuple{Int,Int}}, i::Int)
    @boundscheck checkbounds(g, i)
    return i < length(g) ? g.group_number[1] : g.group_number[2]
end