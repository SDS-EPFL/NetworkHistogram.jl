function create_proposal!(history::MVHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, iteration::Int, proposal.likelihood[1])
    return proposal
end

function make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)
    # copy current in proposal
    deepcopy!(proposal, current)
    # update node labels
    update_labels!(proposal, swap, current)
    # update realized, estimated_theta
    update_observed!(proposal, swap, A)
    # update ll
    updateLL!(proposal, current, swap)
end

function update_labels!(proposal::Assignment, swap::Tuple{Int, Int}, current::Assignment)
    proposal.node_labels[swap[1]] = current.node_labels[swap[2]]
    proposal.node_labels[swap[2]] = current.node_labels[swap[1]]
end

function updateLL!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int})
    # O(G^2) where G is the number of groups
    proposal.likelihood[1] = NetworkHistogram.compute_log_likelihood(proposal)
end

function update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)
    group_updated = [proposal.node_labels[swap[1]], proposal.node_labels[swap[2]]]

    #todo: probability smarter way to update by just summing over the upper triangle
    # of the matrix
    @inbounds for g in 1:length(proposal.group_size)
        for g_prime in group_updated
            proposal.realized[g, g_prime] = sum(A[proposal.node_labels .== g,
                                                  proposal.node_labels .== g_prime])

            # if we look at the connection within the same group, we need to divide by 2
            # to avoid double counting
            if g == g_prime
                proposal.realized[g, g_prime] ÷= 2
            end
            proposal.realized[g_prime, g] = proposal.realized[g, g_prime]
        end
    end

    @. proposal.estimated_theta = proposal.realized / proposal.counts
end
