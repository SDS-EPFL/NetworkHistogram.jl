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

include("swap.jl")
