mutable struct BernoulliSwap{T} <: Swap
    index1::Int
    index2::Int
    old_assignment::BernoulliAssignment{T}
end

function make_swap(assignment::BernoulliAssignment{T}, id::Tuple{Int}) where {T}
    return BernoulliSwap(id[1], id[2], deepcopy(assignment))
end

function make_swap!(swap::BernoulliSwap{T}, assignment::BernoulliAssignment{T},
        id::Tuple{Int}) where {T}
    swap.index1, swap.index2 = id
    swap.old_assignment = deepcopy(assignment)
end

function revert_swap!(assignment::BernoulliAssignment{T}, swap::BernoulliSwap{T}) where {T}
    assignment = deepcopy(swap.old_assignment)
end

function swap!(assignment::BernoulliAssignment{T}, swap::BernoulliSwap{T}) where {T}
    # perform fast update
end
