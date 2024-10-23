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

function deep_copy_matrix_of_vec!(container, source)
    for i in eachindex(container)
        container[i] = copy(source[i])
    end
end

function copy_realized_and_theta!(a,b)
    deep_copy_matrix_of_vec!(a.realized, b.realized)
    deep_copy_matrix_of_vec!(a.estimated_theta, b.estimated_theta)
    a.log_likelihood = b.log_likelihood
    return nothing
end

function make_swap!(
        swap::CategoricalSwap{M, F}, a::CategoricalAssignment{T, M, F, C},
        id::Tuple{Int, Int}) where {T, M, F, C}
    swap.index1, swap.index2 = id
    copy_realized_and_theta!(swap, a.additional_data)
    #copy!.(swap.realized, a.additional_data.realized)
    #copy!.(swap.estimated_theta, a.additional_data.estimated_theta)
    #swap.log_likelihood = a.additional_data.log_likelihood
    #return nothing
end

function revert_swap!(
        a::CategoricalAssignment{T, M, F, C}, swap::CategoricalSwap{M, F}) where {
        T, M, F, C}
    swap_node_labels!(a, swap.index1, swap.index2)
    copy_realized_and_theta!(a.additional_data, swap)
    #copy!.(a.additional_data.realized, swap.realized)
    #copy!.(a.additional_data.estimated_theta, swap.estimated_theta)
    #a.additional_data.log_likelihood = swap.log_likelihood
    #return nothing
end

function apply_swap!(
        a::CategoricalAssignment{T, M, F, C}, swap::CategoricalSwap{M, F}) where {
        T, M, F, C}
    update_observed_and_labels!(a, swap)
    update_ll!(a)
end

function update_ll!(a::CategoricalAssignment)
    a.additional_data.log_likelihood = compute_log_likelihood(
        a.additional_data.estimated_theta, a.additional_data.realized)
    return nothing
end

function fit(
        a::CategoricalAssignment{T, M, F, C}, g::Observations) where {
        T, M, F, C}
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
        a::CategoricalAssignment{T, M, F, C}, swap::CategoricalSwap{M, F}) where {
        T, M, F, C}
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)

    adj_1 = @view a.additional_data.A[:, swap.index1]
    adj_2 = @view a.additional_data.A[:, swap.index2]
    realized_g1 = @view a.additional_data.realized[:, g1]
    realized_g2 = @view a.additional_data.realized[:, g2]

    for i in axes(a.additional_data.A, 1)
        if i == swap.index1 || i == swap.index2
            continue
        end
        obs_1 = adj_1[i]
        obs_2 = adj_2[i]
        group_inter = get_group_of_vertex(a, i)
        if obs_1 != obs_2
            _fast_update!!(a.additional_data.realized, g1, g2, obs_1, obs_2, group_inter)
        end

        # if i == swap.index1 || i == swap.index2 || obs_1 == obs_2
        #     continue
        # else


        #     a_g1_g_inter = a.additional_data.realized[g1, group_inter]
        #     a_g2_g_inter = a.additional_data.realized[g2, group_inter]
        #     a_g_inter_g1 = realized_g1[group_inter]
        #     a_g_inter_g2 = realized_g2[group_inter]

        #     # send from group 1 to group 2
        #     a_g1_g_inter[obs_1] -= 1
        #     a_g_inter_g1[obs_1] = a_g1_g_inter[obs_1]

        #     a_g2_g_inter[obs_1] += 1
        #     a_g_inter_g2[obs_1] = a_g2_g_inter[obs_1]


        #     # send from group 2 to group 1
        #     a_g2_g_inter[obs_2] -= 1
        #     a_g_inter_g2[obs_2] = a_g2_g_inter[obs_2]

        #     a_g1_g_inter[obs_2] += 1
        #     a_g_inter_g1[obs_2] = a_g1_g_inter[obs_2]
        # end
    end

    _fast_div!(a.additional_data.estimated_theta, a.additional_data.realized,
        a.additional_data.counts)

    # swap of the labels should happen after the update of the realized and estimated_theta
    # for the above loop to work correctly
    swap_node_labels!(a, swap.index1, swap.index2)
    return nothing
end


function _fast_update!!(realized, g1, g2, obs_1, obs_2, g_inter)
    realized[g1, g_inter][obs_1] -= 1
    realized[g_inter, g1][obs_1] = realized[g1, g_inter][obs_1]

    realized[g2,g_inter][obs_1] += 1
    realized[g_inter, g2][obs_1] = realized[g2,g_inter][obs_1]

    # send from group 2 to group 1
    realized[g2,g_inter][obs_2] -= 1
    realized[g_inter, g2][obs_2] = realized[g2,g_inter][obs_2]

    realized[g1, g_inter][obs_2] += 1
    realized[g_inter,g1][obs_2] = realized[g1, g_inter][obs_2]
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
