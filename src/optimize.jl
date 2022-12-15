function graphhist(A; h=select_bandwidth(A), maxitr=1000, optimizer = RandomNodeSwap())
    best, current, proposal, history = initialize(A, h, optimizer)

    for i ∈ 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current; )    
        best = update_best!(history, i, current, best)
        if stopping_rule(history, optimizer)
            break
        end
    end

    return GraphHist(best)
end

function create_proposal!(history, i, proposal, current, A, swap_rule)
    swap = select_swap(current, A, swap_rule)
    proposal = make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, i, proposal.likelihood)
    return proposal
end

function update_best!(history, i, current, best)
    if best.likelihood > current.likelihood
        push!(history, :best_likelihood, i, current.likelihood)
        return current
    else
        return best
    end
end

function stopping_rule(history, optimizer)
    # check improvements over last k steps is above threshold
end

function make_proposal!(proposal, current, swap, A)
    # copy current in proposal
    # update based on swap (counts, thetha, ...)
    # update ll
end