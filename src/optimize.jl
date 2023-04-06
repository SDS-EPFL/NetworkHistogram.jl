function checkadjacency(A)
    @assert eltype(A) <: Real
    if !(eltype(A) === Bool)
        @assert all(a ∈ [zero(eltype(A)),one(eltype(A))] for a in A) "All elements of the ajacency matrix should be zero or one."
    end
    @assert issymmetric(A)
    @assert all(A[i,i]==zero(eltype(A)) for i in 1:size(A,1)) "The diagonal of the adjacency matrix should all be zeros."
    return nothing
end

"""
    graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule = RandomNodeSwap(),
    starting_assignment_rule = RandomStart(), accept_rule = Strict(),
    stop_rule = PreviousBestValue(3), record_trace=true)

Compute the graph histogram.

# Arguments
TBW
"""
function graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule::NodeSwapRule = RandomNodeSwap(),
                   starting_assignment_rule::StartingAssignment = RandomStart(), accept_rule::AcceptRule = Strict(),
                   stop_rule::StopRule = PreviousBestValue(3), record_trace = true)

    checkadjacency(A)
    @assert maxitr > 0

    return _graphhist(A, Val{record_trace}(), h = h, maxitr = maxitr, swap_rule = swap_rule,
                      starting_assignment_rule = starting_assignment_rule,
                      accept_rule = accept_rule,
                      stop_rule = stop_rule)
end

"""
    _graphhist(A, record_trace=Val{true}(); h, maxitr, swap_rule, starting_assignment_rule, accept_rule, stop_rule)

Internal version of `graphhist` which is type stable.
"""
function _graphhist(A, record_trace = Val{true}(); h, maxitr, swap_rule,
                    starting_assignment_rule, accept_rule, stop_rule)
    best, current, proposal, history = initialize(A, h, starting_assignment_rule,
                                                  record_trace)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current, accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, stop_rule)
            break
        end
    end

    return graphhist_format_output(best, history)
end

"""
    graphhist_format_output(best, history)

Formates the `graphhist` output depending on the type of `history` requested by the user.
"""
function graphhist_format_output(best, history::TraceHistory)
    return (graphhist = GraphHist(best), trace = history, likelihood = best.likelihood)
end
function graphhist_format_output(best, history::NoTraceHistory)
    return (graphhist = GraphHist(best), likelihood = history.best_likelihood)
end

function update_best!(history::GraphOptimizationHistory, iteration::Int,
                      current::Assignment,
                      best::Assignment)
    if current.likelihood > best.likelihood
        update_best!(history, iteration, current.likelihood)
        deepcopy!(best, current)
    end
    return best
end

"""
    initialize(A, h, starting_assignment_rule, record_trace)

Initialize the memory required for finding optimal graph histogram.
"""
function initialize(A, h, starting_assignment_rule, record_trace)
    node_labels, group_size = initialise_node_labels(A, h, starting_assignment_rule)
    proposal = Assignment(A, node_labels, group_size)
    current = deepcopy(proposal)
    best = deepcopy(proposal)
    history = initialize_history(best, current, proposal, record_trace)
    return best, current, proposal, history
end

function select_bandwith(A)
    error("Automatic bandwidth selection not implemented yet, please specify h manually.")
end
