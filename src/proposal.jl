"""Functions to create and evaluate possible labels update."""

"""
    create_proposal!(history::MVHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)

Create a new proposal by swapping the labels of two nodes. The new assignment is stored in
`proposal`. The swap is selected using the `swap_rule` function. The likelihood of the new
proposal is stored in the history.

!!! warning
    The `proposal` assignment is modified in place to avoid unnecessary memory allocation.
"""
function create_proposal!(history::MVHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, iteration::Int, proposal.likelihood[1])
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
    # update node labels
    update_labels!(proposal, swap, current)
    # update realized, estimated_theta
    update_observed!(proposal, swap, A)
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
    proposal.likelihood[1] = NetworkHistogram.compute_log_likelihood(proposal)
end

"""
    update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)

Update the observed and estimated attributes of the `proposal` assignment based on the
swap of the nodes specified in `swap`.
"""
function update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)
    group_updated = [proposal.node_labels[swap[1]], proposal.node_labels[swap[2]]]

    #todo: probability smarter way to update by just summing over the upper triangle
    # of the matrix
    @inbounds for g in 1:length(proposal.group_size)
        for g_prime in group_updated
            proposal.realized[g, g_prime] = sum(A[proposal.node_labels .== g,
                                                  proposal.node_labels .== g_prime])

            # if we look at the connection within the same group, we need to divide by 2
            # to avoid double counting edges
            if g == g_prime
                proposal.realized[g, g_prime] ÷= 2
            end
            proposal.realized[g_prime, g] = proposal.realized[g, g_prime]
        end
    end

    @. proposal.estimated_theta = proposal.realized / proposal.counts
end
