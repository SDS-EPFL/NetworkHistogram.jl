include("swap_workspace.jl")
include("swap_categorical.jl")
include("config_rules/include.jl")

"""
    GreedyParams

Configuration parameters for the greedy optimization algorithm.

# Fields
- `max_iter::Int`: Maximum number of iterations (default: 100,000)
- `swap_rule::NodeSwapRule`: Rule for selecting which nodes to swap
- `accept_rule::AcceptRule`: Rule for accepting/rejecting proposed swaps
- `stop_rule::StopRule`: Rule for determining when to stop optimization
- `progress_bar::Bool`: Whether to display a progress bar (default: true)

# Examples
```julia
# Use default parameters
params = GreedyParams()

# Custom parameters with stricter stopping
params = GreedyParams(
    1_000_000,                    # max iterations
    RandomNodeSwap(),           # random node selection
    Strict(),                   # only accept improvements
    PreviousBestValue(5000),   # stop after 5000 iterations without improvement
    true                        # show progress bar
)
```

See also: [`NodeSwapRule`](@ref), [`AcceptRule`](@ref), [`StopRule`](@ref)
"""
mutable struct GreedyParams
    max_iter::Int
    swap_rule::NodeSwapRule
    accept_rule::AcceptRule
    stop_rule::StopRule
    progress_bar::Bool
end

"""
    GreedyParams()

Create default greedy optimization parameters.

Defaults:
- max_iter: 100,000
- swap_rule: RandomNodeSwap()
- accept_rule: Strict()
- stop_rule: PreviousBestValue(10,000)
- progress_bar: true
"""
function GreedyParams()
    GreedyParams(
        100_000, RandomNodeSwap(), Strict(), PreviousBestValue(10_000), true)
end

"""
    greedy_optimize(g, initial_labels, params::GreedyParams)

Run greedy optimization to find a good network histogram (block model partition).

# Arguments
- `g`: Tuple of (EdgeList, Dist) containing the network data and distribution type
- `initial_labels`: Initial group assignment for nodes
- `params::GreedyParams`: Optimization parameters

# Returns
- `Assignment`: Optimized assignment of nodes to groups

# Algorithm
The algorithm iteratively:
1. Proposes moving a node to a different group (based on swap_rule)
2. Evaluates the change in log-likelihood
3. Accepts or rejects the move (based on accept_rule)
4. Continues until stopping criterion met (based on stop_rule)
"""
function greedy_optimize(g, initial_labels, params::GreedyParams)
    @debug "making assignment"
    a = Assignment(initial_labels, g...)
    @debug "assignment made, starting greedy search"
    greedy_improve!(a; params = params)
    return a
end

"""
    greedy_improve!(a::Assignment; params = GreedyParams())

Improve an existing assignment through greedy local search.

Modifies the assignment in-place by iteratively proposing and accepting beneficial
node reassignments.

# Arguments
- `a::Assignment`: The assignment to improve (modified in-place)
- `params::GreedyParams`: Optimization parameters

# Note
This function modifies `a` in-place and updates its log-likelihood.
"""
function greedy_improve!(a::Assignment; params = GreedyParams())
    # allocate memory for swap
    swap = make_swap(a, (1, 2))

    # display progress bar
    p = ProgressUnknown(enabled = params.progress_bar,
        showspeed = true, desc = "Greedy search: ")

    for i in 1:(params.max_iter)
        local_search!(a, swap, params)
        next!(p;
            showvalues = [
                ("ll: ", loglikelihood(a)), info_to_print(params.stop_rule)])
        if stopping_rule(a, params.stop_rule)
            if i < 10
                @warn "Greedy search stopped early after $(i) iterations"
            end
            finish!(p)
            break
        end
    end
end

# Internal function for a single local search step
function local_search!(a::Assignment, swap, params::GreedyParams)
    # select two nodes to swap and update data in the swap object
    make_swap!(swap, a, select_indices_swap(a, params.swap_rule))
    # apply swap, test if local improvement and update assignment if needed
    accept_reject_update!(a, swap, params.accept_rule)
end
