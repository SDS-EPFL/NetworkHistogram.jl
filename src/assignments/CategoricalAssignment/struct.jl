mutable struct CategoricalData{M, F}
    counts::Matrix{Int}
    realized::Matrix{MVector{M, Int}}
    estimated_theta::Matrix{MVector{M, F}}
    A::Matrix{Int} # possible use of CategoricalArrays.jl ?
    log_likelihood::F
end

const CategoricalAssignment{T, M, F} = Assignment{T, CategoricalData{M, F}}
const CategoricalInitRule{S, F} = InitRule{S, Val{CategoricalData}}

function CategoricalAssignment(
        g, group_size::GroupSize, node_labels::Vector{Int})
    categorical_data = make_categorical_data(g, node_labels, group_size)
    return Assignment(group_size, node_labels, categorical_data)
end

function make_assignment(g, h, init_rule::CategoricalInitRule)
    group_size,
    node_labels = initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)
    return CategoricalAssignment(g, group_size, node_labels)
end

function make_categorical_data(g, node_labels, group_size)
    number_groups = length(group_size)
    n = length(node_labels)
    A, num_categories = categorical_matrix(g)
    counts = zeros(Int, number_groups, number_groups)
    realized = [MVector{num_categories}(zeros(Int, num_categories))
                for _ in 1:number_groups, _ in 1:number_groups]

    # this is incorrect if the diagonal of the matrix is anything
    # else than 0, and that no "categories" is represented by 0
    @inbounds @simd for k in 1:number_groups
        for l in k:number_groups
            for m in 1:num_categories
                if k == l
                    c = group_size[k] * (group_size[k] - 1) ÷ 2
                    r = sum(A[node_labels .== k, node_labels .== l] .== m) ÷ 2
                else
                    c = group_size[k] * group_size[l]
                    r = sum(A[node_labels .== k, node_labels .== l] .== m)
                end
                realized[k, l][m] = r
                realized[l, k][m] = r
                counts[k, l] = c
                counts[l, k] = c
            end
        end
    end

    estimated_theta = realized ./ counts
    ll = compute_log_likelihood(estimated_theta, realized)
    return CategoricalData(counts, realized, estimated_theta, A, ll)
end

function compute_log_likelihood(
        estimated_theta::AbstractMatrix{MVector{M, T}}, counts::AbstractMatrix{F}) where {
        M, T, F}
    loglik = zero(T)
    number_groups = size(estimated_theta, 1)
    @inbounds for j in 1:number_groups
        @simd for i in j:number_groups
            c = counts[i, j]
            loglik += sum(xlogx.(estimated_theta[i, j]) .* c)
        end
    end
    return loglik
end

# to update, just for test now
function categorical_matrix(A)
    A_inter = A .- minimum(A) .+ 1
    for i in 1:size(A_inter, 1)
        A_inter[i, i] = 0
    end
    return A_inter, maximum(A_inter)
end

function categorical_matrix(g::Observations)
    return categorical_matrix(g.graph)
end

function loglikelihood(a::CategoricalAssignment, g::Observations)
    return a.additional_data.log_likelihood
end

include("swap.jl")
