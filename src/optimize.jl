function graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule = RandomNodeSwap(),
                   starting_assignment_rule = RandomStart(), accept_rule = Strict(),
                   stop_rule = PreviousBestValue(3))
    best, current, proposal, history = initialize(A, h, starting_assignment_rule)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current, accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, stop_rule)
            break
        end
    end

    return GraphHist(best), history
end

function update_best!(history::MVHistory, iteration::Int, current::Assignment,
                      best::Assignment)
    if best.likelihood[1] > current.likelihood[1]
        push!(history, :best_likelihood, iteration::Int, current.likelihood[1])
        return current
    else
        return best
    end
end

function initialize(A, h, starting_assignment_rule)
    node_labels, group_size = initialise_node_labels(A, h, starting_assignment_rule)
    proposal = Assignment(A, node_labels, group_size)
    current = deepcopy(proposal)
    best = deepcopy(proposal)
    history = MVHistory(Dict([
                                 :proposal_likelihood => QHistory(Float64),
                                 :current_likelihood => QHistory(Float64),
                                 :best_likelihood => QHistory(Float64),
                             ]))
    push!(history, :proposal_likelihood, 0, proposal.likelihood[1])
    push!(history, :current_likelihood, 0, current.likelihood[1])
    push!(history, :best_likelihood, 0, best.likelihood[1])
    return best, current, proposal, history
end
