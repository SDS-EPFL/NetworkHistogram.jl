mutable struct SparseSwap{F} <: Swap
    index1::Int
    index2::Int
    realized::Array{Int, 3}
    estimated_theta::Array{F, 3}
    counts::Matrix{Int}
    log_likelihood::F
end

function make_swap(a::SparseAssignment, id)
    return SparseSwap(id[1], id[2], copy(a.additional_data.realized),
        copy(a.additional_data.estimated_theta), copy(a.additional_data.counts),
        a.additional_data.log_likelihood)
end

function copy_addtional!(a, b)
    copy!(a.realized, b.realized)
    copy!(a.estimated_theta, b.estimated_theta)
    copy!(a.counts, b.counts)
    a.log_likelihood = b.log_likelihood
    return nothing
end

function make_swap!(
        swap::SparseSwap{F}, a::SparseAssignment{T, F},
        id) where {T, F}
    swap.index1, swap.index2 = id
    copy_addtional!(swap, a.additional_data)
end

function revert_swap!(
        a::SparseAssignment{T, F}, swap::SparseSwap{F}) where {T, F}
    swap_node_labels!(a, swap.index1, swap.index2)
    copy_addtional!(a.additional_data, swap)
    return nothing
end

function apply_swap!(
        a::SparseAssignment{T, F}, swap::SparseSwap{F}) where {T, F}
    update_observed_and_labels!(a, swap)
    update_ll!(a)
end

function update_ll!(a::SparseAssignment)
    a.additional_data.log_likelihood = compute_log_likelihood_without_0(
        a.additional_data.estimated_theta, a.additional_data.realized, a.additional_data.counts)
    return nothing
end

function update_observed_and_labels!(
        a::SparseAssignment{T, F}, swap::SparseSwap{F}) where {T, F}
    g1 = get_group_of_vertex(a, swap.index1)
    g2 = get_group_of_vertex(a, swap.index2)

    if g1 == g2
        return nothing
    end

    rows = rowvals(a.additional_data.A)
    vals = nonzeros(a.additional_data.A)
    m, n = size(a.additional_data.A)
    for j in [swap.index1, swap.index2]
        a.additional_data.scratch_count .= 0
        a.additional_data.scratch_missing .= 0
        g_from = swap.index1 == j ? g1 : g2
        g_to = swap.index1 == j ? g2 : g1
        for i_index in nzrange(a.additional_data.A, j)
            row = rows[i_index]
            if row == swap.index1 || row == swap.index2
                continue
            end
            val = vals[i_index]
            groupi = get_group_of_vertex(a, row)
            if ismissing(val)
                a.additional_data.scratch_missing[groupi] += 1
            else
                a.additional_data.scratch_count[val, groupi] += 1
            end
        end
        _move_connection!(
            a.additional_data.realized, g_from, g_to, a.additional_data.scratch_count)

        _update_counts!(
            a.additional_data.counts, g_from, g_to, a.additional_data.scratch_missing)
    end

    _fast_div!(a.additional_data.estimated_theta, a.additional_data.realized,
        a.additional_data.counts)

    # swap of the labels should happen after the update of the realized and estimated_theta
    # for the above loop to work correctly
    swap_node_labels!(a, swap.index1, swap.index2)
    return nothing
end

function _update_counts!(counts, g_from, g_to, missing_update)
    for i in axes(counts, 1)
        counts[i, g_to] += missing_update[i]
        counts[i, g_from] -= missing_update[i]
    end
end

function fit(a::SparseAssignment, g::Observations)
    dists = initialize_sbm(a.group_size, ZeroInflatedCategorical(_n_decorations_with_0(a)))
    for group1 in 1:number_groups(a)
        for group2 in 1:number_groups(a)
            theta = a.additional_data.estimated_theta[:, group1, group2]
            dists[group1,
            group2] = ZeroInflatedCategorical(1 - sum(theta), theta)
        end
    end
    return dists
end

function fit(a::SparseAssignment, g::Observations{G, <:DiscretizedDistribution}) where {G}
    dists = initialize_sbm(a.group_size,
        DiscretizedDistribution(
            g.dist_ref.discretizer, ZeroInflatedCategorical(_n_decorations_with_0(a))))
    for group1 in 1:number_groups(a)
        for group2 in 1:number_groups(a)
            theta = a.additional_data.estimated_theta[:, group1, group2]
            p = clamp(1 - sum(theta),0,1)
            dists[group1,
            group2] = DiscretizedDistribution(
                g.dist_ref.discretizer, ZeroInflatedCategorical(p, theta))
        end
    end
    return dists
end
