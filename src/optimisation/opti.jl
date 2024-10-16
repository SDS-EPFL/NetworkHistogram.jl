include("config_rules/swap_rule.jl")
include("config_rules/accept_rule.jl")
include("config_rules/InitRule.jl")
include("config_rules/stop_rule.jl")
include("config_rules/bandwidth_selection_rule.jl")

function greedy_improve!(a::Assignment, G; max_iter::Int = 1000,
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        initialise_rule::InitRule = InitRule(RandomStart(), nothing)
)
end
