function graphhist(A; h = select_bandwidth(A), maxitr = 1000; swap_rule = RandomNodeSwap(), starting_assignment_rule, accept_rule, stop_rule)
    best, current, proposal, history = initialize(A, h; starting_assignment_rule=starting_assignment_rule)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current; accept_rule = accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history; stop_rule=stop_rule)
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