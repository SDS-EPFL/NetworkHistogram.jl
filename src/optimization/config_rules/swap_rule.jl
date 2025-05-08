abstract type NodeSwapRule end

struct RandomNodeSwap <: NodeSwapRule end
struct RandomGroupSwap <: NodeSwapRule end
"""
    select_indices_swap(node_assignment::Assignment, ::NodeSwapRule)

Selects two nodes to swap based on the `NodeSwapRule`, the adjacency matrix `A` and the
current assignment `node_assignment`.

# Implemented rules
- `RandomNodeSwap()`: Select two nodes at random.
- `RandomGroupSwap()`: Select two nodes from two different groups at random.
"""
select_swap

function select_indices_swap(assignment::Assignment, ::RandomNodeSwap)
    return StatsBase.sample(1:number_nodes(assignment), 2; replace = false)
end

function select_indices_swap(assignment::Assignment, ::RandomGroupSwap)
    groups = StatsBase.sample(
        1:number_groups(assignment), 2; replace = false)
    index1 = rand(findall(x -> x == groups[1], assignment.node_labels))
    index2 = rand(findall(x -> x == groups[2], assignment.node_labels))
    return (index1, index2)
end
