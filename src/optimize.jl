"""
    graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule = RandomNodeSwap(),
    starting_assignment_rule = RandomStart(), accept_rule = Strict(),
    stop_rule = PreviousBestValue(3), record_trace=true)

Computes the graph histogram approximation.

# Arguments
- `A`: adjacency matrix of a simple graph

- `h`: bandwidth of the graph histogram (number of nodes in a group or percentage (in [0,1]) of
    nodes in a group)

- `record_trace` (optional): whether to record the trace of the optimization process and return
    it as part of the output. Default is `true`.

# Returns
named tuple with the following fields:
- `graphhist`: the graph histogram approximation
- `trace`: the trace of the optimization process (if `record_trace` is `true`)
- `likelihood`: the loglikelihood of the graph histogram approximation

# Examples
```julia
julia> A = [0 0 1 0 1 0 1 1 0 1
             0 0 1 1 1 1 1 1 0 0
             1 1 0 1 0 0 0 0 1 0
             0 1 1 0 1 0 1 0 0 0
             1 1 0 1 0 0 1 0 0 1
             0 1 0 0 0 0 0 1 0 0
             1 1 0 1 1 0 0 1 0 1
             1 1 0 0 0 1 1 0 0 1
             0 0 1 0 0 0 0 0 0 1
             1 0 0 0 1 0 1 1 1 0]
julia> out = graphhist(A);
julia> graphist_approx = out.graphhist
...
julia> trace = out.trace
NetworkHistogram.TraceHistory{...}
  :best_likelihood => 1 elements {Int64,Float64}
  :proposal_likelihood => 5 elements {Int64,Float64}
  :current_likelihood => 5 elements {Int64,Float64})
julia> loglikelihood = out.likelihood
-22.337057781338277
```
"""
function graphhist(A; h = select_bandwidth(A), maxitr = 10000,
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        starting_assignment_rule::StartingAssignment = EigenStart(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(100), record_trace = true,
        show_progress = true)
    checkadjacency(A)
    @assert maxitr > 0
    A = drop_disconnected_components(A)

    return _graphhist(A, Val{record_trace}(), h = h, maxitr = maxitr, swap_rule = swap_rule,
        starting_assignment_rule = starting_assignment_rule,
        accept_rule = accept_rule,
        stop_rule = stop_rule,
        show_progress = show_progress)
end

"""
    _graphhist(A, record_trace=Val{true}(); h, maxitr, swap_rule, starting_assignment_rule, accept_rule, stop_rule)

Internal version of `graphhist` which is type stable.
"""
function _graphhist(A, record_trace = Val{true}(); h, maxitr, swap_rule,
        starting_assignment_rule, accept_rule, stop_rule, show_progress = true)
    best, current, proposal, history, A = initialize(A, h, starting_assignment_rule,
        record_trace)
    p = Progress(length(maxitr); enabled = show_progress)
    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current, accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, stop_rule)
            break
        end
        next!(p)
    end

    if show_progress
        finish!(p)
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
    node_labels, group_size = initialize_node_labels(A, h, starting_assignment_rule)
    proposal = Assignment(A, node_labels, group_size)
    current = deepcopy(proposal)
    best = deepcopy(proposal)
    history = initialize_history(best, current, proposal, record_trace)
    return best, current, proposal, history, update_adj(A)
end


function update_best!(history::GraphOptimizationHistory, iteration::Int,
        current::Assignments,
        best::Assignments)
    if current.likelihood > best.likelihood
        update_best!(history, iteration, current.likelihood)
        deepcopy!(best, current)
    end
    return best
end



function graphhist(A, d; h = select_bandwidth(A, d), maxitr = 10000,
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        starting_assignment_rule::StartingAssignment = RandomStart(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(100), record_trace = true,
        show_progress = true)
    @assert maxitr > 0
    A = drop_disconnected_components(A)
    best, current, proposal, history, A = initialize(A, d, h, starting_assignment_rule,
        Val{record_trace}())
    p = Progress(length(maxitr); enabled = show_progress)
    for i in 1:maxitr
        swap = select_swap(current, A, swap_rule)
        proposal = update(current, A, swap)
        current = accept_reject_update!(history, i, proposal, current, accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, stop_rule)
            break
        end
        next!(p)
    end
    if show_progress
        finish!(p)
    end
    return best, history
end
