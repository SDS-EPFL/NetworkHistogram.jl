abstract type GraphOptimizationHistory end
struct TraceHistory{M <: MVHistory} <: GraphOptimizationHistory
    history::M
end
mutable struct NoTraceHistory <: GraphOptimizationHistory
    current_iteration::Int
    best_iteration::Int
    best_likelihood::Float64
end

"""
    initialize_history(best, current, proposal, ::Val{true})

    initialize the history when `record_trace=true` is passed to `graphhist`.
"""
function initialize_history(best, current, proposal, ::Val{true})
    history = MVHistory(Dict([
        :proposal_likelihood => QHistory(Float64),
        :current_likelihood => QHistory(Float64),
        :best_likelihood => QHistory(Float64),
    ]))
    push!(history, :proposal_likelihood, 0, proposal.likelihood)
    push!(history, :current_likelihood, 0, current.likelihood)
    push!(history, :best_likelihood, 0, best.likelihood)
    return TraceHistory(history)
end

"""
    initialize_history(best, current, proposal, ::Val{false})

initialize the history when `record_trace=false` is passed to `graphhist`.
"""
function initialize_history(best, current, proposal, ::Val{false})
    return NoTraceHistory(0, 0, best.likelihood)
end

"""
    get_currentitr(history::GraphOptimizationHistory)

Return the current iteration of the optimization from the history.
"""
get_currentitr(history::TraceHistory) = last(history.history, :current_likelihood)
get_currentitr(history::NoTraceHistory) = history.current_iteration

"""
    get_bestitr(history::GraphOptimizationHistory)

Return the best iteration of the optimization from the history.
"""
get_bestitr(history::TraceHistory) = last(history.history, :best_likelihood)
get_bestitr(history::NoTraceHistory) = history.best_iteration

"""
    update_current!(history::GraphOptimizationHistory, iteration, likelihood)

Updates the current value and iteration in history.
"""
function update_current!(history::TraceHistory, iteration, likelihood)
    push!(history.history, :current_likelihood, iteration, likelihood)
end
function update_current!(history::NoTraceHistory, iteration, likelihood)
    history.current_iteration = iteration
end

"""
    update_best!(history::GraphOptimizationHistory, iteration, likelihood)

Updates the best value and iteration in history.
"""
function update_best!(history::TraceHistory, iteration, likelihood)
    push!(history.history, :best_likelihood, iteration, likelihood)
end
function update_best!(history::NoTraceHistory, iteration, likelihood)
    history.best_iteration = iteration
    history.best_likelihood = likelihood
end

"""
    update_previous!(history::GraphOptimizationHistory, iteration, likelihood)

Updates the previous value and iteration in history.

Note this does not apply is `history` is a `NoTraceHistory`, so nothign happens.
"""
function update_proposal!(history::TraceHistory, iteration, likelihood)
    push!(history.history, :proposal_likelihood, iteration, likelihood)
end
update_proposal!(history::NoTraceHistory, iteration, likelihood) = nothing
