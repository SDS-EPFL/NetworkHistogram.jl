mutable struct BernoulliData{F}
    counts::Matrix{Int}
    realized::Matrix{Int}
    estimated_theta::Matrix{F}
    A::BitMatrix  # possible improvement by using an adjacency list  Graphs.SimpleGraphs.adj(G)
    log_likelihood::F
end

const BernoulliAssignment{T, F} = Assignment{T, BernoulliData{F}}
const BernoulliInitRule{S, F} = InitRule{S, Val{BernoulliData}}

function BernoulliAssignment(
        g, group_size::GroupSize, node_labels::Vector{Int})
    bernoulli_data = make_bernoulli_data(g, node_labels, group_size)
    return Assignment(group_size, node_labels, bernoulli_data)
end

function make_assignment(g, h, init_rule::BernoulliInitRule)
    group_size,
    node_labels = initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)
    return BernoulliAssignment(g, group_size, node_labels)
end

# might be worth using graph accessors instead of the adjacency matrix ?
function make_bernoulli_data(g, node_labels, group_size)
    number_groups = length(group_size)
    n = length(node_labels)
    counts = zeros(Int, number_groups, number_groups)
    realized = zeros(Int, number_groups, number_groups)
    A = convert_bitmatrix(g)

    # below needs to be abstracted: not sure how diagonal is handled if nonzero
    # addtioally, we should be able to deal with missing values !
    # This concerns the counts matrix above as well
    @inbounds @simd for k in 1:number_groups
        for l in k:number_groups
            realized[k, l] = sum(A[node_labels .== k, node_labels .== l])
            realized[l, k] = realized[k, l]
            counts[k, l] = group_size[k] * group_size[l]
            counts[l, k] = counts[k, l]
        end
    end

    @inbounds @simd for k in 1:number_groups
        counts[k, k] = group_size[k] * (group_size[k] - 1) ÷ 2
        realized[k, k] = sum(A[node_labels .== k, node_labels .== k]) ÷ 2
    end

    estimated_theta = realized ./ counts
    ll = compute_log_likelihood(estimated_theta, counts)
    return BernoulliData(counts, realized, estimated_theta, A, ll)
end

function convert_bitmatrix(g::Observations{<:AbstractGraph, D}) where {D}
    A = collect(adjacency_matrix(g.graph))
    return convert(BitMatrix, collect(adjacency_matrix(g.graph)))
end

function convert_bitmatrix(g::Observations{<:AbstractMatrix, D}) where {D}
    return convert(BitMatrix, g.graph)
end

function compute_log_likelihood(estimated_theta::AbstractMatrix{F},
        counts::AbstractMatrix{T}) where {F <: Real, T <: Real}
    number_groups = size(estimated_theta, 1)
    loglik = zero(eltype(estimated_theta))
    @inbounds for j in 1:number_groups
        @simd for i in j:number_groups
            θ = estimated_theta[i, j]
            loglik += (xlogx(θ) + xlogx(1 - θ)) * counts[i, j]
        end
    end
    return loglik
end

function log_likelihood(assignment::BernoulliAssignment)
    return assignment.additional_data.log_likelihood
end

log_likelihood(a::BernoulliAssignment, g::Observations) = log_likelihood(a)

function force_recompute_ll(a::BernoulliAssignment, g::Observations)
    a_simple = Assignment(a.group_size, a.node_labels)
    return log_likelihood(a_simple, g)
end

function get_ordered_adjacency_matrix(a::BernoulliAssignment)
    perm = sortperm(a.node_labels)
    return a.additional_data.A[perm, perm]
end

include("swap.jl")
