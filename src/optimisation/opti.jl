include("config_rules/include.jl")

function optimize(g, h = select_bandwidth(g);
        max_iter::Int = 1000,
        initialise_rule::InitRule = InitRule(RandomStart(), nothing),
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(10),
        progress_bar::Bool = false
)
    a = make_assignment(g, h, initialise_rule)
    initialise_stop_rule!(stop_rule, a, g)
    greedy_improve!(a, g; max_iter, swap_rule, accept_rule, stop_rule, progress_bar)
    return a
end

function greedy_improve!(a::Assignment, g; max_iter::Int = 1000,
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(10),
        progress_bar::Bool = false,
)
    # swap memory allocation
    swap = make_swap(a, (1, 1))
    p = Progress(max_iter; enabled = progress_bar)
    # perform local search until the stopping rule is met
    for i in 1:max_iter
        local_search!(a, g, swap, swap_rule = swap_rule, accept_rule = accept_rule)
        next!(p)
        if stopping_rule(a, g, stop_rule)
            finish!(p)
            break
        end
    end
end

# perform local search by trying a swap and accepting it if it improves the likelihood
function local_search!(a::Assignment, g, swap::Swap = make_swap(a, (1, 1));
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict()
)
    # select two nodes to swap and build the swap object
    make_swap!(swap, a, select_swap(a, swap_rule))
    # perform the swap and accept it if it improves the likelihood
    accept_reject_update!(a, swap, g, accept_rule)
end
