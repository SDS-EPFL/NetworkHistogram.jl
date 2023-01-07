function graphhist(A::Matrix{Int}; h = select_bandwidth(A), maxitr = 1000, optimizer = RandomNodeSwap())
    best, current, proposal, history = initialize(A, h, optimizer)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A)
        current = accept_reject_update!(history, i, proposal, current;)
        best = update_best!(history, i, current, best)
        if stopping_rule(proposal, current, best;)
            break
        end
    end

    return GraphHist(best)
end

function create_proposal!(history::MultivalueHistory,
                          iteration::Int,
                          proposal::Assignment,
                          current::Assignment,
                          A::Matrix{Int})
    swap = select_swap(current, A)
    make_proposal!(proposal, current, swap, A)
    push!(history, :proposal_likelihood, iteration, proposal.likelihood)
    return proposal
end

function update_best!(history::MultivalueHistory,
                      iteration::Int,
                      current::Assignment,
                      best::Assignment)
    if best.likelihood > current.likelihood
        push!(history, :best_likelihood, iteration, current.likelihood)
        return current
    else
        return best
    end
end

function stopping_rule(history, optimizer)
    # check improvements over last k steps is above threshold
end
