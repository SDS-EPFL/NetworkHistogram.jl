abstract type StopRule end

function initialise_stop_rule!(stop_rule::StopRule, a, g)
end

# default score is the log likelihood
function score(a::Assignment, g::Observations)
    return log_likelihood(a, g) / binomial(number_nodes(a), 2)
end

mutable struct PreviousBestValue{T} <: StopRule
    k::Int
    past_values::CircularDeque{T}
    function PreviousBestValue(
            k::Int, x::T = -Inf) where {T <: Real}
        @assert k > 0
        # queue stores the best values and at most k subsequent values
        queue = CircularDeque{T}(k + 1)
        push!(queue, x)
        new{T}(k, queue)
    end
end

function initialise_stop_rule!(stop_rule::PreviousBestValue, a, g)
    score_value = score(a, g)
    empty!(stop_rule.past_values)
    push!(stop_rule.past_values, score_value)
end

"""
    stopping_rule(assignment::Assignment,g, stop_rule::StopRule)

Returns a Bool with true if we should stop the optimization based on the `stop_rule`.

# Implemented rules
- `PreviousBestValue(k)`: Stop if the current iteration is `k` iterations away from the
  iteration with the best value.
"""
stopping_rule

function stopping_rule(assignment::Assignment, g, stop_rule::PreviousBestValue)
    score_value = score(assignment, g)
    if isempty(stop_rule.past_values)
        push!(stop_rule.past_values, score_value)
        return false
    elseif score_value > first(stop_rule.past_values)
        empty!(stop_rule.past_values)
        push!(stop_rule.past_values, score_value)
        return false
    else
        return length(stop_rule.past_values) == capacity(stop_rule.past_values)
    end
end
