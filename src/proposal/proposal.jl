function make_proposal!(proposal, current, swap, A)
    # copy current in proposal
    # update based on swap (counts, thetha, ...)
    # update ll
end

function create_proposal!(history, i, proposal, current, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    proposal = make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, i, proposal.likelihood)
    return proposal
end
