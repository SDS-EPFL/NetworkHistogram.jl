function make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)
    # copy current in proposal
    # update based on swap (counts, thetha, ...)
    # update ll
end

function create_proposal!(history::MVHistory, iteration::Int, proposal::Assignment,
                          current::Assignment, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    proposal = make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, iteration::Int, proposal.likelihood)
    return proposal
end
