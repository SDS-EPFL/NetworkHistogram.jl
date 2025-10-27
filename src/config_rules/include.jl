include("swap_rule.jl")
include("stop_rule.jl")

struct GreedyParams{N <: NodeSwapRule, S <: StopRule}
    max_iter::Int
    node_swap_rule::N
    stop_rule::S
end
