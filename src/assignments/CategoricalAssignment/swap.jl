mutable struct CategoricalSwap{M, F} <: Swap
    index1::Int
    index2::Int
    realized::Matrix{MVector{M, Int}}
    estimated_theta::Matrix{MVector{M, F}}
    log_likelihood::F
end

function make_swap(a::CategoricalAssignment, id::Tuple{Int, Int})
    return CategoricalSwap(id[1], id[2], copy(a.additional_data.realized),
        copy(a.additional_data.estimated_theta),
        a.additional_data.log_likelihood)
end

function make_swap!(
        swap::CategoricalSwap{M, F}, a::CategoricalAssignment{T, M, F},
        id::Tuple{Int, Int}) where {T, M, F}
    swap.index1, swap.index2 = id
    copy!(swap.realized, a.additional_data.realized)
    copy!(swap.estimated_theta, a.additional_data.estimated_theta)
    swap.log_likelihood = a.additional_data.log_likelihood
end

function revert_swap!(
        a::CategoricalAssignment{T, M, F}, swap::CategoricalSwap{M, F}) where {
        T, M, F}
    swap_node_labels!(a, swap.index1, swap.index2)
    copy!(a.additional_data.realized, swap.realized)
    copy!(a.additional_data.estimated_theta, swap.estimated_theta)
    a.additional_data.log_likelihood = swap.log_likelihood
end

function apply_swap!(
        a::CategoricalAssignment{T, M, F}, swap::CategoricalSwap{M, F}) where {
        T, M, F}
    update_observed_and_labels!(a, swap)
    update_ll!(a)
end

function update_ll!(a::CategoricalAssignment)
    a.additional_data.log_likelihood = compute_log_likelihood(
        a.additional_data.estimated_theta, a.additional_data.counts)
    return nothing
end

function fit_sbm(
        a::CategoricalAssignment{T, M, F}, g::Observations) where {T, M, F}
    dists = initialize_sbm(a.group_size, Categorical(ones(M) / M))
    for group1 in 1:number_groups(a)
        for group2 in 1:number_groups(a)
            dists[group1,
            group2] = Categorical(a.additional_data.estimated_theta[
                group1, group2])
        end
    end
    return dists
end

function update_observed_and_labels!(
        a::CategoricalAssignment{T, M, F}, swap::CategoricalSwap{M, F}) where {
        T, M, F}
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)

    adj_1 = @view a.additional_data.A[:, swap.index1]
    adj_2 = @view a.additional_data.A[:, swap.index2]
    realized_g1 = @view a.additional_data.realized[:, g1]
    realized_g2 = @view a.additional_data.realized[:, g2]

    @inbounds @fastmath for i in axes(a.additional_data.A, 1)
        index_1 = adj_1[i]
        index_2 = adj_2[i]
        if i == swap.index1 || i == swap.index2 || index_1 == index_2

        else
            group_inter = get_group_of_vertex(a, i)

            a_g1_g_inter = a.additional_data.realized[g1, group_inter]
            a_g2_g_inter = a.additional_data.realized[g2, group_inter]
            a_g_inter_g1 = realized_g1[group_inter]
            a_g_inter_g2 = realized_g2[group_inter]

            # send from group 1 to group 2
            a_g1_g_inter[index_1] -= 1
            a_g_inter_g1[index_1] -= 1

            a_g2_g_inter[index_2] += 1
            a_g_inter_g2[index_2] += 1

            # send from group 2 to group 1
            a_g2_g_inter[index_2] -= 1
            a_g_inter_g2[index_2] -= 1

            a_g1_g_inter[index_1] += 1
            a_g_inter_g1[index_1] += 1
        end
    end

    _fast_div!(a.additional_data.estimated_theta, a.additional_data.realized,
        a.additional_data.counts)

    # swap of the labels should happen after the update of the realized and estimated_theta
    # for the above loop to work correctly
    swap_node_labels!(a, swap.index1, swap.index2)
    return nothing
end

function _fast_div!(theta, realized, counts)
    for j in axes(theta, 2)
        for i in axes(theta, 1)
            t = theta[i, j]
            for k in axes(t, 1)
                theta[i, j][k] = realized[i, j][k] / counts[i, j]
            end
        end
    end
end
