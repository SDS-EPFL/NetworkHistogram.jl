abstract type Swap end

mutable struct DefaultSwap <: Swap
    index1::Int
    index2::Int
end

function make_swap(::Assignment{T, Nothing}, id::Tuple{Int, Int}) where {T}
    return DefaultSwap(id[1], id[2])
end

function make_swap!(swap::DefaultSwap, a::Assignment, id::Tuple{Int, Int})
    swap.index1, swap.index2 = id
end

function apply_swap!(a::Assignment, s::DefaultSwap)
    a.node_labels[s.index1], a.node_labels[s.index2] = a.node_labels[s.index2],
    a.node_labels[s.index1]
end

revert_swap!(assignment::Assignment, swap::DefaultSwap) = apply_swap!(assignment, swap)
