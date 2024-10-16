include("config_rules/swap_rule.jl")
include("config_rules/accept_rule.jl")
include("config_rules/InitRule.jl")
include("config_rules/stop_rule.jl")
include("config_rules/bandwidth_selection_rule.jl")

function optimize(G, h; initialise_rule::InitRule = InitRule(RandomStart(), nothing))
    a = make_assignment(G, h, initialise_rule)
    greedy_improve!(a, G)
    return a
end

function optimize(G, h = select_bandwidth(G);
        max_iter::Int = 1000,
        initialise_rule::InitRule = InitRule(RandomStart(), nothing),
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = MaxIter(1000)
)
    a = make_assignment(G, h, initialise_rule)
    for i in 1:max_iter
        local_search!(a, G, swap_rule = swap_rule, accept_rule = accept_rule)
        if stop_rule(a, G, stop_rule)
            break
        end
    end
    greedy_improve!(a, G; max_iter, swap_rule, accept_rule, stop_rule)
    return a
end

function greedy_improve!(a::Assignment, G, max_iter::Int = 1000;
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = MaxIter(1000)
)
end

# perform local search by trying a swap and accepting it if it improves the likelihood
function local_search!(a::Assignment, G;
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict()
)
    swap = select_swap(a, swap_rule)
    if accept_swap(a, swap, G, accept_rule)
        apply_swap!(a, swap)
    end
end
