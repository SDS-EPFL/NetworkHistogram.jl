include("swap_rule.jl")
include("stop_rule.jl")

abstract type ParamsType end

@kwdef struct GreedyParams{N <: NodeSwapRule, S <: StopRule} <: ParamsType
    max_iter::Int = 1_000_000
    stalled_iters::Int = 5_000
    node_swap_rule::N = RandomGroupSwap()
    stop_rule::S = PreviousBestValue(stalled_iters, Inf, :min)
    display_progress::Bool = true
    progress_freq::Int = 10_000
    warm_start::Bool = false
end

function reset!(params::GreedyParams)
    reset!(params.stop_rule)
    return params
end
