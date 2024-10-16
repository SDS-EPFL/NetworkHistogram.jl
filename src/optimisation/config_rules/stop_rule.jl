abstract type StopRule end
mutable struct PreviousBestValue{T} <: StopRule
    k::Int
    past_values::Queue{T}
    function PreviousBestValue(
            k::Int, x::T=0.0) where {T <: Real}
        @assert k > 0
        new{T}(k, Queue{T}())
    end
end

"""
    stopping_rule(assignment::Assignment,G, stop_rule::StopRule)

Returns a Bool with true if we should stop the optimization based on the `stop_rule`.

# Implemented rules
- `PreviousBestValue(k)`: Stop if the current iteration is `k` iterations away from the
  iteration with the best value.
"""
stopping_rule

function stopping_rule(assignment::Assignment, G, stop_rule::PreviousBestValue)
    log_likelihood = log_likelihood(assignment, G)
    if length(stop_rule.past_values) == 0
        push!(stop_rule.past_values, log_likelihood)
        return false
    elseif log_likelihood > first(stop_rule.past_values)
        empty!(stop_rule.past_values)
        push!(stop_rule.past_values, log_likelihood)
        return false
    else
        push!(stop_rule.past_values, log_likelihood)
        return length(stop_rule.past_values) == stop_rule.k + 1 #always keep the best value
    end
end
