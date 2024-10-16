abstract type StopRule end

function initialise_stop_rule!(stop_rule::StopRule, a, g)
end

# default score is the log likelihood
function score(a::Assignment, g::Observations)
    return log_likelihood(a, g) / binomial(number_nodes(a), 2)
end

mutable struct PreviousBestValue{T} <: StopRule
    k::Int
    previous_best_value::T
    iterations_since_best::Int
    function PreviousBestValue(
            k::Int, x::T = -Inf) where {T <: Real}
        @assert k > 0
        # queue stores the best values and at most k subsequent values
        new{T}(k, x, 0)
    end
end

function initialise_stop_rule!(stop_rule::PreviousBestValue, a, g)
    score_value = score(a, g)
    stop_rule.previous_best_value = score_value
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
    if score_value > stop_rule.previous_best_value
        stop_rule.previous_best_value = score_value
        stop_rule.iterations_since_best = 0
    else
        stop_rule.iterations_since_best += 1
    end
    return stop_rule.iterations_since_best >= stop_rule.k
end
