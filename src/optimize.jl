function graphhist(A; h = select_bandwidth(A), maxitr = 1000, optimizer = RandomNodeSwap())
    best, current, proposal, history = initialize(A, h, optimizer)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current;)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, optimizer)
            break
        end
    end

    return GraphHist(best)
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

function accept_reject_update!(history, i, proposal, current;)
    # fill in

    push!(history, :likelihood, i, current.likelihood)
    return current
end
