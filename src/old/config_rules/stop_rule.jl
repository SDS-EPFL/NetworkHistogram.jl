abstract type StopRule end
struct PreviousBestValue <: StopRule
    k::Int
    function PreviousBestValue(k::Int)
        @assert k > 0
        new(k)
    end
end

"""
    stopping_rule(history, stop_rule::StopRule)

Returns a Bool with true if we should stop the optimization based on the `stop_rule`.

# Implemented rules
- `PreviousBestValue(k)`: Stop if the current iteration is `k` iterations away from the
  iteration with the best value.
"""
stopping_rule

function stopping_rule(history::GraphOptimizationHistory, stop_rule::PreviousBestValue)
    current_itr = get_currentitr(history)
    prev_best_itr = get_bestitr(history)
    return current_itr[1] - prev_best_itr[1] > stop_rule.k
end
