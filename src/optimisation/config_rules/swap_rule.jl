abstract type NodeSwapRule end

struct RandomNodeSwap <: NodeSwapRule end

"""
    select_swap(node_assignment::Assignment, ::NodeSwapRule)

Selects two nodes to swap based on the `NodeSwapRule`, the adjacency matrix `A` and the
current assignment `node_assignment`.

# Implemented rules
- `RandomNodeSwap()`: Select two nodes at random.
"""
select_swap

function select_swap(assignment::Assignment, ::RandomNodeSwap)
    groups = StatsBase.sample(
        1:number_groups(assignment), 2; replace = false)
    index1 = rand(get_vertex_in_group(assignment, groups[1]))
    index2 = rand(get_vertex_in_group(assignment, groups[2]))
    return (index1, index2)
end
