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

function SparseAssignment( g::Observations{G,D}, group_size::GroupSize, node_labels::Vector{Int}) where {G,D}
    A = issparse(g.graph) ? g.graph : sparse(g.graph)
    num_levels = length(unique(A)) -1
    sparse_data = SparseData(A, size(group_size, 1), num_levels, group_size, node_labels)
    return Assignment(group_size, node_labels, sparse_data)
end


function SparseData(A::SparseMatrixCSC{T, Int}, k::Int,
        level_count::Int, group_size, node_labels) where {T}
    n = size(A, 1)
    data = SparseData(zeros(Int, k, k), zeros(Int, level_count, k, k),
        zeros(Float64, level_count, k, k), dropzeros!(A), zeros(Int, level_count, k), zeros(Int, k), 0.0)
    _count_possible_occurences!(data, group_size)
    _count_occurences!(data, node_labels)
    _fast_div!(data.estimated_theta, data.realized, data.counts)
    println(data.estimated_theta)
    println(data.realized)
    println(data.counts)
    data.log_likelihood = compute_log_likelihood_without_0(data.estimated_theta, data.realized, data.counts)
    return data
end


function _count_possible_occurences!(data, group_size)
    k = size(group_size, 1)
    for j in 1:k
        data.counts[j, j] = group_size[j] * (group_size[j] - 1) ÷ 2
        for i in j+1:k
            data.counts[i, j] = group_size[i] * group_size[j]
            data.counts[j, i] = group_size[i] * group_size[j]
        end
    end
end

function _count_occurences!(data, node_labels)
    rows = rowvals(data.A)
    vals = nonzeros(data.A)
    m, n = size(data.A)
    for j in 1:n
        groupj = node_labels[j]
        for i in nzrange(data.A, j)
            row = rows[i]
            val = vals[i]
            groupi = node_labels[row]
            if ismissing(val)
                data.counts[groupj, groupj] -= 1
                if groupj != groupj
                    data.counts[groupj, groupj] -= 1
                end
            else
                data.realized[val, groupi, groupj] += 1
                if groupi != groupj
                    data.realized[val, groupj, groupi] += 1
                end
            end
        end
    end
end



function compute_log_likelihood_without_0(
        estimated_theta::Array{T, 3}, realized::Array{F, 3}, counts) where {
        T, F}
    loglik = zero(T)
    number_groups = size(estimated_theta, 2)
    number_decorations = size(estimated_theta, 1)
    @inbounds for j in 1:number_groups
        for i in j:number_groups
            prob_absent = one(T)
            total_decorations = counts[i, j]
            for m in 1:number_decorations
                if realized[m, i, j] != 0
                    prob_absent -= estimated_theta[m, i, j]
                    total_decorations -= realized[m, i, j]
                    loglik += realized[m, i, j] * log(estimated_theta[m, i, j])
                end
            end
            println(total_decorations, prob_absent)
            loglik += total_decorations * log(prob_absent)
        end
    end
    return loglik
end


include("swap.jl")
