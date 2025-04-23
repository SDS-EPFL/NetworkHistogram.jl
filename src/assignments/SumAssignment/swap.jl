mutable struct SumSwap{F} <: Swap
    index1::Int
    index2::Int
    θ::Dict{Tuple{Int, Int}, F}
    counts::Dict{Tuple{Int, Int}, Int}
    log_likelihood_per_group::Dict{Tuple{Int, Int}, Float64}
    log_likelihood::Float64
end

function make_swap(a::SumAssignment, id)
    return SumSwap(id[1], id[2], deepcopy(a.additional_data.θ),
        deepcopy(a.additional_data.counts), deepcopy(a.additional_data.log_likelihood_per_group), a.additional_data.log_likelihood)
end


function make_swap!(swap::SumSwap{F}, a::SumAssignment{T, F}, id) where {T, F}
    swap.index1, swap.index2 = id
    swap.θ = deepcopy(a.additional_data.θ)
    swap.counts = deepcopy(a.additional_data.counts)
    swap.log_likelihood_per_group = deepcopy(a.additional_data.log_likelihood_per_group)
    swap.log_likelihood = a.additional_data.log_likelihood
end

function revert_swap!(
        a::SumAssignment{T, F}, swap::SumSwap{F}) where {T, F}
    swap_node_labels!(a, swap.index1, swap.index2)
    a.additional_data.θ = deepcopy(swap.θ)
    a.additional_data.counts = deepcopy(swap.counts)
    a.additional_data.log_likelihood_per_group = deepcopy(swap.log_likelihood_per_group)
    a.additional_data.log_likelihood = swap.log_likelihood
end

function apply_swap!(
        a::SumAssignment{T, F}, swap::SumSwap{F}) where {T, F}
    λ = a.additional_data.λ
    rows = rowvals(λ)
    vals = nonzeros(λ)
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)
    if g1 == g2
        return nothing
    end

    for v in rows[nzrange(λ, swap.index1)]
        key_old_groups = minmax(g1, a.node_labels[v])
        key_new_groups = minmax(g2, a.node_labels[v])
        c_og = a.counts[key_old_groups]
        c_ng = a.counts[key_new_groups]
        param = vals[i]
        a.θ[key_old_groups] = (a.θ[key_old_groups]*c_og - param)/(c_og - 1)
        a.θ[key_new_groups] = (a.θ[key_new_groups]*c_ng + param)/(c_ng + 1)
        a.counts[key_old_groups] -= 1
        a.counts[key_new_groups] += 1
    end

    for v in rows[nzrange(λ, swap.index2)]
        key_old_groups = minmax(g2, a.node_labels[v])
        key_new_groups = minmax(g1, a.node_labels[v])
        c_og = a.counts[key_old_groups]
        c_ng = a.counts[key_new_groups]
        param = vals[i]
        a.θ[key_old_groups] = (a.θ[key_old_groups]*c_og - param)/(c_og - 1)
        a.θ[key_new_groups] = (a.θ[key_new_groups]*c_ng + param)/(c_ng + 1)
        a.counts[key_old_groups] -= 1
        a.counts[key_new_groups] += 1
    end

    swap_node_labels!(a, swap.index1, swap.index2)
    fast_update_ll!(a, swap)
end

function fast_update_ll(a::SumAssignment, swap::SumSwap)
    k = size(a.group_size, 1)
    for i in 1:k
        for j in i:k
            index_group = (i, j)
            if swap.θ[index_group] != a.θ[index_group]
                _update_ll_one_group!(a, index_group)
            end
        end
    end
    a.additional_data.log_likelihood = sum(values(a.additional_data.log_likelihood_per_group))
end

function _update_ll_one_group!(a::SumAssignment, group)
    k = size(a.group_size, 1)
    nodes_1 = findall(x -> x == group[1], a.node_labels)
    nodes_2 = findall(x -> x == group[2], a.node_labels)
    ll = 0.0
    rows = rowvals(a.additional_data.λ)
    vals = nonzeros(a.additional_data.λ)
    for i in nodes_1
        for u in nodes_1
            for v in rows[nzrange(a.additional_data.λ, u)]
                if v ∈ nodes_2
                    ll += loglikelihood(a.θ[group], a.additional_data.A[u, v])
                end
            end
        end
    end
    a.log_likelihood_per_group[group] = ll
    return nothing
end
