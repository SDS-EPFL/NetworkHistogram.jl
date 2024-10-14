abstract type NodeSwapRule end

struct RandomNodeSwap <: NodeSwapRule end

"""
    select_swap(node_assignment::Assignment, A, ::NodeSwapRule)

Selects two nodes to swap based on the `NodeSwapRule`, the adjacency matrix `A` and the
current assignment `node_assignment`.

# Implemented rules
- `RandomNodeSwap()`: Select two nodes at random.
"""
select_swap

function select_swap(node_assignment::Assignment, A, ::RandomNodeSwap)
    index1 = rand(1:size(A, 1))
    label1 = node_assignment.node_labels[index1]
    index2 = index1
    for _ in 1:10
        index2 = rand(1:size(A, 1))
        if node_assignment.node_labels[index2] != label1
            break
        end
    end
    return (index1, index2)
end
