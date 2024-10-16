abstract type Swap end

mutable struct DefaultSwap <: Swap
    index1::Int
    index2::Int
end

function make_swap(::Assignment, id::Tuple{Int, Int})
    return DefaultSwap(id[1], id[2])
end

function make_swap!(swap::DefaultSwap, a::Assignment, id::Tuple{Int, Int})
    swap.index1, swap.index2 = id
end

apply_swap!(a::Assignment, s::DefaultSwap) = swap_node_labels!(a, s.index1, s.index2)

revert_swap!(assignment::Assignment, swap::DefaultSwap) = apply_swap!(assignment, swap)

function swap_node_labels!(a::Assignment, i, j)
    a.node_labels[i], a.node_labels[j] = a.node_labels[j],a.node_labels[i]
end
