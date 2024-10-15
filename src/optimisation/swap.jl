abstract type Swap end

mutable struct DefaultSwap <: Swap
    index1::Int
    index2::Int
end

make_swap(::Assignment{T,Nothing}, id::Tuple{Int}) where {T} = DefaultSwap(id[1], id[2])
function make_swap!(swap::DefaultSwap, assignment::Assignment, id::Tuple{Int})
    swap.index1, swap.index2 = id
end

function swap!(assignment::Assignment, swap::DefaultSwap)
    assignment.node_labels[swap.index1], assignment.node_labels[swap.index2] = assignment.node_labels[swap.index2],
    assignment.node_labels[swap.index1]
end

revert_swap!(assignment::Assignment, swap::DefaultSwap) = swap!(assignment, swap)
