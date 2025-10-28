abstract type StopRule end

function info_to_print(::StopRule)
    return nothing
end

mutable struct PreviousBestValue{T, S} <: StopRule
    const k::Int
    previous_best_value::T
    iterations_since_best::Int
end

function PreviousBestValue(k::Int, x::T = -Inf, best = :max) where {T <: Real}
    @argcheck k > 0
    PreviousBestValue{T, Val(best)}(k, x, 0)
end

const PreviousMaxValue{T} = PreviousBestValue{T, Val(:max)}
const PreviousMinValue{T} = PreviousBestValue{T, Val(:min)}

function reset!(stop_rule::PreviousBestValue{T}, loss_value::T) where {T}
    stop_rule.previous_best_value = loss_value
    stop_rule.iterations_since_best = 0
end

reset!(stop_rule::PreviousMaxValue) = reset!(stop_rule, -Inf)
reset!(stop_rule::PreviousMinValue) = reset!(stop_rule, Inf)

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
        reset!(stop_rule, loss)
    else
        stop_rule.iterations_since_best += 1
    end
    return stop_rule.iterations_since_best >= stop_rule.k
end

function info_to_print(stop_rule::PreviousBestValue)
    ("stalled iter: ", stop_rule.iterations_since_best)
end
