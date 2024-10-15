abstract type StopRule end
mutable struct PreviousBestValue{T} <: StopRule
    k::Int
    past_values::Queue{T}
    function PreviousBestValue(k::Int)
        @assert k > 0
        new{T}(k, Queue{T}(k))
    end
end

"""
    stopping_rule(assignment::Assignment, stop_rule::StopRule)

Returns a Bool with true if we should stop the optimization based on the `stop_rule`.

# Implemented rules
- `PreviousBestValue(k)`: Stop if the current iteration is `k` iterations away from the
  iteration with the best value.
"""
stopping_rule

function stopping_rule(assignment::Assignment, stop_rule::PreviousBestValue)
    loglikelihood = loglikelihood(assignment)
    if length(stop_rule.past_values) == 0
        push!(stop_rule.past_values, loglikelihood)
        return false
    elseif loglikelihood > first(stop_rule.past_values)
        empty!(stop_rule.past_values)
        push!(stop_rule.past_values, loglikelihood)
        return false
    else
        push!(stop_rule.past_values, loglikelihood)
        return length(stop_rule.past_values) == stop_rule.k + 1 #always keep the best value
    end
end
