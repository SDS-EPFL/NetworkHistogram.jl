abstract type StopRule end

function info_to_print(::StopRule)
    return nothing
end

function initialise_stop_rule!(stop_rule::StopRule, a, g)
end

function score(a::Assignment)
    return loglikelihood(a)
end

mutable struct PreviousBestValue{T, S} <: StopRule
    k::Int
    previous_best_value::T
    iterations_since_best::Int
end

function PreviousBestValue(k::Int, x::T = -Inf, best = :max) where {T <: Real}
    @argcheck k > 0
    PreviousBestValue{T, Val(best)}(k, x, 0)
end

const PreviousMaxValue{T} = PreviousBestValue{T, Val(:max)}
const PreviousMinValue{T} = PreviousBestValue{T, Val(:min)}

function initialise_stop_rule!(stop_rule::PreviousBestValue, a)
    score_value = score(a)
    stop_rule.previous_best_value = score_value
end

function compare_to_best(current, past, ::PreviousMaxValue)
    return current > past
end

function compare_to_best(current, past, ::PreviousMinValue)
    return current < past
end

"""
    stopping_rule(assignment::Assignment,g, stop_rule::StopRule)

Returns a Bool with true if we should stop the optimization based on the `stop_rule`.

# Implemented rules
- `PreviousBestValue(k)`: Stop if the current iteration is `k` iterations away from the
  iteration with the best value.
"""
stopping_rule

function stopping_rule(loss::T, stop_rule::PreviousBestValue{T}) where {T <: Real}
    if compare_to_best(loss, stop_rule.previous_best_value, stop_rule)
        stop_rule.previous_best_value = loss
        stop_rule.iterations_since_best = 0
    else
        stop_rule.iterations_since_best += 1
    end
    return stop_rule.iterations_since_best >= stop_rule.k
end

stopping_rule(a, stop_rule::StopRule) = stopping_rule(score(a), stop_rule)

function info_to_print(stop_rule::PreviousBestValue)
    ("stalled iter: ", stop_rule.iterations_since_best)
end
