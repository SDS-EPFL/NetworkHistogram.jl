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

function update_best!(history::MVHistory, iteration::Int, current::Assignment,
                      best::Assignment)
    if best.likelihood > current.likelihood
        push!(history, :best_likelihood, iteration::Int, current.likelihood)
        return current
    else
        return best
    end
end

function stopping_rule(history::MVHistory, optimizer)
    # check improvements over last k steps is above threshold
end

function accept_reject_update!(history::MVHistory, iteration::Int, proposal::Assignment,
                               current::Assignment;)
    # fill in

    push!(history, :likelihood, iteration, current.likelihood)
    return current
end
