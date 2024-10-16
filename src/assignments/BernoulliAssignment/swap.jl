mutable struct BernoulliSwap{F} <: Swap
    index1::Int
    index2::Int
    realized::Matrix{Int}
    estimated_theta::Matrix{F}
    log_likelihood::F
    node_labels::Vector{Int}
end

function make_swap(assignment::BernoulliAssignment{T, F}, id::Tuple{Int, Int}) where {T, F}
    return BernoulliSwap(id[1], id[2], copy(assignment.additional_data.realized),
        copy(assignment.additional_data.estimated_theta),
        assignment.additional_data.log_likelihood, copy(assignment.node_labels))
    # realized = copy(assignment.additional_data.realized)
    # estimated_theta = copy(assignment.additional_data.estimated_theta)
    # log_likelihood = assignment.additional_data.log_likelihood
    # return BernoulliSwap(id[1], id[2], realized, estimated_theta, log_likelihood)
end

function make_swap!(swap::BernoulliSwap{F}, assignment::BernoulliAssignment{T, F},
        id::Tuple{Int, Int}) where {T, F}
    swap.index1, swap.index2 = id
    copy!(swap.realized, assignment.additional_data.realized)
    copy!(swap.estimated_theta, assignment.additional_data.estimated_theta)
    #copy!(swap.node_labels, assignment.node_labels)
    swap.log_likelihood = assignment.additional_data.log_likelihood
end

function revert_swap!(
        assignment::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    swap_node_labels!(assignment, swap.index1, swap.index2)
    copy!(assignment.additional_data.realized, swap.realized)
    copy!(assignment.additional_data.estimated_theta, swap.estimated_theta)
    #copy!(assignment.node_labels, swap.node_labels)
    assignment.additional_data.log_likelihood = swap.log_likelihood
end

function apply_swap!(
        assignment::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    # swap of the labels should happen after the update of the realized and estimated_theta
    update_observed!(assignment, swap)
    swap_node_labels!(assignment, swap.index1, swap.index2)
    update_ll!(assignment)
end

function update_observed!(a::BernoulliAssignment{T, F}, swap::BernoulliSwap{F}) where {T, F}
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)

    for i in axes(a.additional_data.A, 2)
        if i == swap.index1 || i == swap.index2 ||
           a.additional_data.A[swap.index1, i] == a.additional_data.A[swap.index2, i]
            continue
        end
        group_inter = get_group_of_vertex(a, i)
        if a.additional_data.A[swap.index1, i]
            a.additional_data.realized[g1, group_inter] -= 1
            a.additional_data.realized[group_inter, g1] = a.additional_data.realized[
                g1, group_inter]

            a.additional_data.realized[g2, group_inter] += 1
            a.additional_data.realized[group_inter, g2] = a.additional_data.realized[
                g2, group_inter]
        end
        if a.additional_data.A[swap.index2, i]
            a.additional_data.realized[g2, group_inter] -= 1
            a.additional_data.realized[group_inter, g2] = a.additional_data.realized[
                g2, group_inter]

            a.additional_data.realized[g1, group_inter] += 1
            a.additional_data.realized[group_inter, g1] = a.additional_data.realized[
                g1, group_inter]
        end
    end

    @. a.additional_data.estimated_theta = a.additional_data.realized /
                                           a.additional_data.counts
    return nothing
end

function update_ll!(a::BernoulliAssignment)
    a.additional_data.log_likelihood = compute_log_likelihood(
        a.additional_data.estimated_theta, a.additional_data.counts)
    return nothing
end


function fit(a::BernoulliAssignment, g::Observations)
    println("Fitting BernoulliAssignment")
    dists = initialize_sbm(a.group_size, Bernoulli(0.5))
    for group1 in 1:number_groups(a)
        for group2 in 1:number_groups(a)
            dists[group1, group2] = Bernoulli(a.additional_data.estimated_theta[group1, group2])
        end
    end
    return dists
end
