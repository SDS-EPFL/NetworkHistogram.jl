struct BernoulliData{T}
    counts::Matrix{Int}
    realized::Matrix{Int}
    estimated_theta::Matrix{T}
    A::BitMatrix
end

const BernoulliAssignment{T} = Assignment{T, BernoulliData}
const BernoulliInitRule{S} = InitRule{S, Val{BernoulliData}}


# is this type stable? should this be BernoulliAssignment{T,F}? see line 8 above
function BernoulliAssignment(
        G, node_labels::Vector{Int}, group_size::GroupSize{T}) where {T}
    k = length(group_size)
    return BernoulliAssignment{T}(group_size, node_labels,
        BernoulliData(zeros(Int, k, k), zeros(Int, k, k), zeros(T, k, k), BitMatrix(G)))
end

function make_assignment(G, h, init_rule::BernoulliInitRule{S}) where {S}
    return BernoulliAssignment(initialize_node_labels(
        G, h, init_rule.starting_assignment_rule)...)
end

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



function test_init_rule(rule::InitRule{S, Val{T}}) where {S, T <: BernoulliData}
    println("Method with BernoulliData")
    println("Starting assignment rule: ", rule.starting_assignment_rule)
end


function test_init_rule_with_alias(rule::BernoulliInitRule{S}) where {S}
    println("Method with BernoulliData")
    println("Starting assignment rule: ", rule.starting_assignment_rule)
end
