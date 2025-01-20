mutable struct SparseData{F, C}
    counts::Matrix{Int}
    realized::Array{Int, 3}
    estimated_theta::Array{F, 3}
    A::SparseMatrixCSC{C, Int}
    scratch_count::Matrix{Int}
    scratch_missing::Vector{Int}
    log_likelihood::F
end

const SparseAssignment{T, F, C} = Assignment{
    T, SparseData{F, C}}
const SparseInitRule{S, F} = InitRule{S, Val{SparseData}}

function SparseAssignment(
        g::Observations{G, D}, group_size::GroupSize, node_labels::Vector{Int}) where {
        G, D}
    A = issparse(g.graph) ? g.graph : sparse(g.graph)
    num_levels = ncategories(g.dist_ref)
    sparse_data = SparseData(
        A, size(group_size, 1), num_levels, group_size, node_labels)
    return Assignment(group_size, node_labels, sparse_data)
end

function make_assignment(g, h, init_rule::SparseInitRule)
    group_size,
    node_labels = initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)
    return SparseAssignment(g, group_size, node_labels)
end

function SparseData(A::SparseMatrixCSC{T, Int}, k::Int,
        level_count::Int, group_size, node_labels) where {T}
    n = size(A, 1)
    data = SparseData(zeros(Int, k, k), zeros(Int, level_count, k, k),
        zeros(Float64, level_count, k, k), dropzeros(A), zeros(
            Int, level_count, k), zeros(
            Int, k), 0.0)
    _count_possible_occurences!(data, group_size)
    _count_occurences!(data, node_labels)
    _fast_div!(data.estimated_theta, data.realized, data.counts)
    data.log_likelihood = compute_log_likelihood_without_0(
        data.estimated_theta, data.realized, data.counts)
    return data
end

function _count_possible_occurences!(data, group_size)
    k = size(group_size, 1)
    for j in 1:k
        data.counts[j, j] = group_size[j] * (group_size[j] - 1) ÷ 2
        for i in (j + 1):k
            data.counts[i, j] = group_size[i] * group_size[j]
            data.counts[j, i] = group_size[i] * group_size[j]
        end
    end
end

function _count_occurences!(data, node_labels)
    m, n = size(data.A)
    for k in 1:length(unique(node_labels))
        for l in k:length(unique(node_labels))
            node_group_k = findall(x -> x == k, node_labels)
            node_group_l = findall(x -> x == l, node_labels)
            if k != l
                counts = StatsBase.countmap(data.A[i, j] for i in node_group_k
                for j in node_group_l if i != j)
            else
                counts = StatsBase.countmap(data.A[i, j] for i in node_group_k
                for j in node_group_l if i < j)
            end
            for m in 1:size(data.realized, 1)
                data.realized[m, k, l] = get(counts, m, 0)
                data.realized[m, l, k] = get(counts, m, 0)
            end
            total_witouth_missing = sum(values(counts)) -
                                    get(counts, missing, 0)
            data.counts[k, l] = total_witouth_missing
            data.counts[l, k] = total_witouth_missing
        end
    end
end

function compute_log_likelihood_without_0(
        estimated_theta::Array{T, 3}, realized::Array{F, 3}, counts) where {
        T, F}
    loglik = zero(T)
    number_groups = size(estimated_theta, 2)
    number_decorations = size(estimated_theta, 1)
    for j in 1:number_groups
        for i in j:number_groups
            total_decorations = counts[i, j]
            if total_decorations < sum(realized[:, i, j])
                total_decorations = sum(realized[:, i, j])
            end
            loglik -= xlogx(total_decorations)
            for m in 1:number_decorations
                loglik += xlogx(realized[m, i, j])
                total_decorations -= realized[m, i, j]
            end
            loglik += xlogx(total_decorations)
        end
    end
    return loglik
end

function _n_decorations_with_0(a::SparseAssignment)
    return size(a.additional_data.estimated_theta, 1) + 1
end

function _n_decorations_not_0(a::SparseAssignment)
    return size(a.additional_data.estimated_theta, 1)
end

function loglikelihood(assignment::SparseAssignment, g::Observations)
    return assignment.additional_data.log_likelihood
end

include("swap.jl")
