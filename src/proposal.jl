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
    #todo: check how to bound away from 0 and 1
    group1 = proposal.node_labels[swap[1]]
    group2 = proposal.node_labels[swap[2]]

    for g in 1:length(proposal.group_size)
        proposal.likelihood[1] += proposal.counts[g, group1] *
                                  (log(proposal.realized[g, group1]) -
                                   log(current.realized[g, group1]))
        proposal.likelihood[1] += proposal.counts[g, group2] *
                                  (log(proposal.realized[g, group2]) -
                                   log(current.realized[g, group2]))
        proposal.likelihood[1] += proposal.realized[g, group1] - current.realized[g, group1]
        proposal.likelihood[1] += proposal.realized[g, group2] - current.realized[g, group2]
    end
end

function update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)
    group1 = proposal.node_labels[swap[1]]
    group2 = proposal.node_labels[swap[2]]

    #todo: probability smarter way to update by just summing over the upper triangle
    # of the matrix
    for g in 1:length(proposal.group_size)
        proposal.realized[g, group1] = sum(A[proposal.node_labels .== g,
                                             proposal.node_labels .== group1]) ÷ 2
        proposal.realized[g, group2] = sum(A[proposal.node_labels .== g,
                                             proposal.node_labels .== group2]) ÷ 2
        proposal.realized[group1, g] = proposal.realized[g, group1]
        proposal.realized[group2, g] = proposal.realized[g, group2]
    end

    proposal.estimated_theta .= proposal.realized ./ proposal.counts
end
