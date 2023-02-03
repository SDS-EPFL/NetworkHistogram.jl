abstract type StopRule end
struct PreviousBestValue
    k::Int
    function PreviousBestValue(k::Int)
        @assert k > 0
        new(k)
    end
end

function stopping_rule(history::MVHistory; stop_rule::PreviousBestValue)
    current_itr = last(history, :current_likelihood)
    prev_best_itr = last(history, :best_likelihood)
    return current_itr - prev_best_itr > stop_rule.k
end