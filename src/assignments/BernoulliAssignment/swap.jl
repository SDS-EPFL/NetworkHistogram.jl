mutable struct BernoulliSwap{F} <: Swap
    index1::Int
    index2::Int
    realized::Matrix{Int}
    estimated_theta::Matrix{F}
    log_likelihood::F
    node_labels::Vector{Int}
end

function make_swap(
        a::BernoulliAssignment{T, F}, id::Tuple{Int, Int}) where {T, F}
    return BernoulliSwap(id[1], id[2], copy(a.additional_data.realized),
        copy(a.additional_data.estimated_theta),
        a.additional_data.log_likelihood, copy(a.node_labels))
end

function make_swap!(swap::BernoulliSwap{F}, a::BernoulliAssignment{T, F},
        id::Tuple{Int, Int}) where {T, F}
    swap.index1, swap.index2 = id
    copy!(swap.realized, a.additional_data.realized)
    copy!(swap.estimated_theta, a.additional_data.estimated_theta)
    swap.log_likelihood = a.additional_data.log_likelihood
end

function revert_swap!(
        a::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    swap_node_labels!(a, swap.index1, swap.index2)
    copy!(a.additional_data.realized, swap.realized)
    copy!(a.additional_data.estimated_theta, swap.estimated_theta)
    a.additional_data.log_likelihood = swap.log_likelihood
end

function apply_swap!(
        a::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    update_observed_and_labels!(a, swap)
    update_ll!(a)
end

function update_observed_and_labels!(
        a::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)

    for i in axes(a.additional_data.A, 2)
        if i == swap.index1 || i == swap.index2 ||
           a.additional_data.A[swap.index1, i] ==
           a.additional_data.A[swap.index2, i]
            continue
        end
        group_inter = get_group_of_vertex(a, i)
        if a.additional_data.A[swap.index1, i]
            a.additional_data.realized[g1, group_inter] -= 1
            a.additional_data.realized[
                group_inter, g1] = a.additional_data.realized[
                g1, group_inter]

            a.additional_data.realized[g2, group_inter] += 1
            a.additional_data.realized[
                group_inter, g2] = a.additional_data.realized[
                g2, group_inter]
        end
        if a.additional_data.A[swap.index2, i]
            a.additional_data.realized[g2, group_inter] -= 1
            a.additional_data.realized[
                group_inter, g2] = a.additional_data.realized[
                g2, group_inter]

            a.additional_data.realized[g1, group_inter] += 1
            a.additional_data.realized[
                group_inter, g1] = a.additional_data.realized[
                g1, group_inter]
        end
    end

    @. a.additional_data.estimated_theta = a.additional_data.realized /
                                           a.additional_data.counts

    # swap of the labels should happen after the update of the realized and estimated_theta
    # for the above loop to work correctly
    swap_node_labels!(a, swap.index1, swap.index2)
    return nothing
end

function update_ll!(a::BernoulliAssignment)
    a.additional_data.log_likelihood = compute_log_likelihood(
        a.additional_data.estimated_theta, a.additional_data.counts)
    return nothing
end

function fit_sbm(a::BernoulliAssignment, g::Observations)
    dists = initialize_sbm(a.group_size, Bernoulli(0.5))
    for group1 in 1:number_groups(a)
        for group2 in 1:number_groups(a)
            dists[group1,
                group2] = Bernoulli(a.additional_data.estimated_theta[
                group1, group2])
        end
    end
    return dists
end
