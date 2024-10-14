"""Functions to create and evaluate possible labels update."""

"""
    create_proposal!(history::GraphOptimizationHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)

Create a new proposal by swapping the labels of two nodes. The new assignment is stored in
`proposal`. The swap is selected using the `swap_rule` function. The likelihood of the new
proposal is stored in the history.

!!! warning
    The `proposal` assignment is modified in place to avoid unnecessary memory allocation.
"""
function create_proposal!(history::GraphOptimizationHistory, iteration::Int,
        proposal::Assignment,
        current::Assignment, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    make_proposal!(proposal, current, swap, A)
    update_proposal!(history, iteration, proposal.likelihood)
    return proposal
end

"""
    make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)

From the current assignment, create a new assignment by swapping the labels of the nodes
specified in `swap`. The new assignment is stored in `proposal`.
"""
function make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)
    # copy current in proposal
    deepcopy!(proposal, current)
    # update realized, estimated_theta
    update_observed!(proposal, swap, A)
    # update node labels (has to happen after!!!)
    update_labels!(proposal, swap, current)
    # update ll
    updateLL!(proposal)
end

"""
    update_labels!(proposal::Assignment, swap::Tuple{Int, Int}, current::Assignment)

Update the labels of the nodes specified in `swap` in the `proposal` assignment.
"""
function update_labels!(proposal::Assignment, swap::Tuple{Int, Int}, current::Assignment)
    proposal.node_labels[swap[1]] = current.node_labels[swap[2]]
    proposal.node_labels[swap[2]] = current.node_labels[swap[1]]
end

"""
    updateLL!(proposal::Assignment)

Update the likelihood of the `proposal` assignment based on its observed and estimated
attributes.
"""
function updateLL!(proposal::Assignment)
    # O(G^2) where G is the number of groups
    proposal.likelihood = NetworkHistogram.compute_log_likelihood(proposal)
end

"""
    update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)

Update the observed and estimated attributes of the `proposal` assignment based on the
swap of the nodes specified in `swap`.

NOTE labels of the nodes before the swap
"""

function update_observed!(proposal::Assignment{T, 2}, swap::Tuple{Int, Int}, A) where {T}
    group_node_1 = proposal.node_labels[swap[1]]
    group_node_2 = proposal.node_labels[swap[2]]

    for i in axes(A, 1)
        if i == swap[1] || i == swap[2] || A[swap[1], i] == A[swap[2], i]
            continue
        end
        group_i = proposal.node_labels[i]
        if A[i, swap[1]] == 1
            proposal.realized[group_node_1, group_i] -= 1
            proposal.realized[group_i, group_node_1] = proposal.realized[group_node_1,
                group_i]

            proposal.realized[group_node_2, group_i] += 1
            proposal.realized[group_i, group_node_2] = proposal.realized[group_node_2,
                group_i]
        end
        if A[i, swap[2]] == 1
            proposal.realized[group_node_2, group_i] -= 1
            proposal.realized[group_i, group_node_2] = proposal.realized[group_node_2,
                group_i]

            proposal.realized[group_node_1, group_i] += 1
            proposal.realized[group_i, group_node_1] = proposal.realized[group_node_1,
                group_i]
        end
    end

    @. proposal.estimated_theta = proposal.realized / proposal.counts

    return nothing
end

function update_observed!(proposal::Assignment{T, 3}, swap::Tuple{Int, Int}, A) where {T}
    group_node_1 = proposal.node_labels[swap[1]]
    group_node_2 = proposal.node_labels[swap[2]]
    if group_node_1 == group_node_2
        return nothing
    end

    for i in axes(A, 1)
        if i == swap[1] || i == swap[2] || A[swap[1], i] == A[swap[2], i]
            continue
        end
        group_i = proposal.node_labels[i]

        proposal.realized[group_node_1, group_i, A[i, swap[1]]] -= 1
        proposal.realized[group_i, group_node_1, A[i, swap[1]]] = proposal.realized[group_node_1,
            group_i, A[i, swap[1]]]
        proposal.realized[group_node_2, group_i, A[i, swap[1]]] += 1
        proposal.realized[group_i, group_node_2, A[i, swap[1]]] = proposal.realized[group_node_2,
            group_i, A[i, swap[1]]]

        proposal.realized[group_node_1, group_i, A[i, swap[2]]] += 1
        proposal.realized[group_i, group_node_1, A[i, swap[2]]] = proposal.realized[group_node_1,
            group_i, A[i, swap[2]]]
        proposal.realized[group_node_2, group_i, A[i, swap[2]]] -= 1
        proposal.realized[group_i, group_node_2, A[i, swap[2]]] = proposal.realized[group_node_2,
            group_i, A[i, swap[2]]]
    end

    @. proposal.estimated_theta = proposal.realized / proposal.counts

    return nothing
end
