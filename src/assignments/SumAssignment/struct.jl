
# type F needs to be a vector field!

struct SumData{F, C}
    λ::SparseMatrixCSC{F, Int}
    θ::Dict{Tuple{Int, Int}, F}
    A::SparseMatrixCSC{C, Int}
    counts::Dict{Tuple{Int, Int}, Int}
    log_likelihood_per_group::Dict{Tuple{Int, Int}, Float64}
    log_likelihood::Float64
end

const SumAssignment{T, F, C} = Assignment{T, SumData{F, C}}
const SumInitRule{S} = InitRule{S, Val{SumData}}

function SumAssignment(
        A::SparseMatrixCSC{C, Int},
        λ::SparseMatrixCSC{F, Int}, group_size::GroupSize, node_labels::Vector{Int}) where {
        F, C}
    k = size(group_size, 1)
    θ = Dict{Tuple{Int, Int}, F}()
    counts = Dict{Tuple{Int, Int}, Int}()

    rows = rowvals(λ)
    vals = nonzeros(λ)
    m, n = size(λ)
    for u in 1:n
        for v in rows[nzrange(λ, u)]
            if u >= v
                break # check that this isn't a mistake trying to be fast
                continue
            end
            key_groups = minmax(node_labels[u], node_labels[v])
            param = vals[i]
            if haskey(θ, key_groups)
                θ[key_groups] += param
            else
                θ[key_groups] = param
            end
            if haskey(counts, key_groups)
                counts[key_groups] += 1
            else
                counts[key_groups] = 1
            end
        end
    end
    for i in 1:k
        for j in i:k
            θ[minmax(i, j)] ./= counts[minmax(i, j)]
        end
    end
    for i in 1:k
        counts[(i, i)] ./= 2
    end
    ll_sum = 0.0
    ll = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:k
        for j in i:k
            ll[(i, j)] = 0.0
        end
    end
    for i in 1:n
        for v in nzrange(λ, j)
            u = rows[i]
            if u >= v
                continue
            end
            key_groups = minmax(node_labels[u], node_labels[v])
            ll[minmax(
                node_labels[u], node_labels[v])] += loglikelihood(θ[key_groups], A[u, v])
        end
    end
    ll_sum = sum(values(ll))
    return Assignment(group_size, node_labels, SumData(λ, θ, A, counts, ll, ll_sum))
end

function loglikelihood(assignment::SumAssignment, g::Observations)
    return sum(values(assignment.additional_data.log_likelihood))
end

include("swap.jl")
