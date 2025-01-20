include("config_rules/include.jl")

"""
    estimate_graphon(graph, h; iterations, initialise_rule, swap_rule, accept_rule, stop_rule, progress_bar)

Estimate the graphon for the given graph.

# Arguments
- `graph`: The input graph.
- `h`: Number of nodes per block.
- `iterations`: Maximum number of iterations.
- `initialise_rule::InitRule`: Rule for initializing the assignment.
- `swap_rule::NodeSwapRule`: Rule for swapping nodes.
- `accept_rule::AcceptRule`: Rule for accepting swaps.
- `stop_rule::StopRule`: Rule for stopping the iterations.
- `progress_bar::Bool`: Whether to show a progress bar.

# Returns
- `a`: The assignment of nodes to blocks.
"""
function estimate_graphon(
        graph, h = select_number_node_per_block(graph, EstimatedDegrees());
        iterations::Int = 10_000,
        initialise_rule::InitRule = InitRule(SpectralStart(), nothing),
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(1000),
        progress_bar::Bool = false
)
    a = make_assignment(graph, h, initialise_rule)
    initialise_stop_rule!(stop_rule, a, graph)
    greedy_improve!(
        a, graph; iterations, swap_rule, accept_rule, stop_rule, progress_bar)
    return a
end

"""
    greedy_improve!(a::Assignment, graph; iterations, swap_rule, accept_rule, stop_rule, progress_bar)

Perform greedy improvement on the assignment.

# Arguments
- `a::Assignment`: The assignment of nodes to blocks.
- `graph`: The input graph.
- `iterations`: Maximum number of iterations.
- `swap_rule::NodeSwapRule`: Rule for swapping nodes.
- `accept_rule::AcceptRule`: Rule for accepting swaps.
- `stop_rule::StopRule`: Rule for stopping the iterations.
- `progress_bar::Bool`: Whether to show a progress bar.
"""
function greedy_improve!(a::Assignment, graph; iterations::Int = 10_000,
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict(),
        stop_rule::StopRule = PreviousBestValue(1000),
        progress_bar::Bool = false
)
    # swap memory allocation
    swap = make_swap(a, (1, 1))
    p = ProgressUnknown(
        enabled = progress_bar, showspeed = true, desc = "Greedy search: ")
    # perform local search until the stopping rule is met
    for i in 1:iterations
        local_search!(
            a, graph, swap, swap_rule = swap_rule, accept_rule = accept_rule)
        next!(p)
        if stopping_rule(a, graph, stop_rule)
            finish!(p)
            break
        end
    end
end

"""
    local_search!(a::Assignment, graph, swap; swap_rule, accept_rule)

Perform local search by trying a swap and accepting it if it improves the likelihood.

# Arguments
- `a::Assignment`: The assignment of nodes to blocks.
- `graph`: The input graph.
- `swap`: The swap object.
- `swap_rule::NodeSwapRule`: Rule for swapping nodes.
- `accept_rule::AcceptRule`: Rule for accepting swaps.
"""
function local_search!(
        a::Assignment, graph, swap::Swap = make_swap(a, (1, 1));
        swap_rule::NodeSwapRule = RandomNodeSwap(),
        accept_rule::AcceptRule = Strict()
)
    # select two nodes to swap and build the swap object
    make_swap!(swap, a, select_swap(a, swap_rule))
    # perform the swap and accept it if it improves the likelihood
    accept_reject_update!(a, swap, graph, accept_rule)
end
