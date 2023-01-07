abstract type NodeSwapRule end

struct RandomNodeSwap <: NodeSwapRule end

function select_swap(node_assignment::Assignment, A, method::RandomNodeSwap)::Tuple{Int, Int}
    group_index = sort(sample(1:(node_assignment.number_groups), 2, replace = false)) # make non-allocating version of this
    node_index1 = sample(1:node_assignment.group_size[1])
    n2 = group_index == node_assignment.number_groups ? node_assignment.group_size[2] :
         node_assignment.group_size[1]
    node_index2 = sample(1:n2)
    return (node_index1, node_index2)
end
