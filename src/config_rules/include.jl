include("swap_rule.jl")
include("stop_rule.jl")

abstract type ParamsType end

@kwdef struct GreedyParams{N <: NodeSwapRule, S <: StopRule} <: ParamsType
    max_iter::Int = 1_000_000
    stalled_iters::Int = 5_000
    node_swap_rule::N = RandomGroupSwap()
    stop_rule::S = PreviousBestValue(stalled_iters, Inf, :min)
end

function reset!(params::GreedyParams)
    params.stop_rule = PreviousBestValue(
        params.stalled_iters, Inf, :min)
    return params
end
