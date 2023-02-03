function create_proposal!(history::MVHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    proposal = make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, iteration::Int, proposal.likelihood[1])
    return proposal
end

function make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)
    # copy current in proposal
    deepcopy!(proposal, current)
    # update node labels
    proposal.node_labels[swap[1]] = current.node_labels[swap[2]]
    proposal.node_labels[swap[2]] = current.node_labels[swap[1]]

    # update counts, realized, estimated_theta
    #update_observed!(proposal, swap, A)
    # update ll
    #proposal = Assignment(proposal,updateLL(proposal, current,A))
    # for now just create new assignment from scratch
    proposal = Assignment(A, proposal.node_labels, current.group_size)
end

function updateLL(proposal::Assignment, current::Assignment, A)
end

function update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)
end
