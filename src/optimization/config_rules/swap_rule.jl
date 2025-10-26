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

function select_indices_swap(node_labels::AbstractVector{Int}, ::RandomNodeSwap)
    return Tuple(StatsBase.samplepair(1:length(node_labels)))
end

function select_indices_swap(node_labels::AbstractVector{Int}, ::RandomGroupSwap,
        k::Int = length(unique(node_labels)))
    groups = StatsBase.sample(1:k, 2; replace = false)
    index1 = rand(findall(x -> x == groups[1], node_labels))
    index2 = rand(findall(x -> x == groups[2], node_labels))
    return index1, index2
end

function select_indices_swap(a::Assignment, rule::NodeSwapRule)
    select_indices_swap(a.node_labels, rule)
end

function select_indices_swap(assignment::Assignment, rule::RandomGroupSwap)
    return select_indices_swap(assignment.node_labels, rule, number_groups(assignment))
end
