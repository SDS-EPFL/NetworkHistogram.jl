include("swap_workspace.jl")
include("swap_categorical.jl")
include("config_rules/include.jl")

mutable struct GreedyParams
    max_iter::Int
    swap_rule::NodeSwapRule
    accept_rule::AcceptRule
    stop_rule::StopRule
    progress_bar::Bool
end

function GreedyParams()
    GreedyParams(
        100_000, RandomNodeSwap(), Strict(), PreviousBestValue(10_000), true)
end

function greedy_optimize(g, initial_labels, params::GreedyParams)
    @debug "making assignment"
    a = Assignment(initial_labels, g...)
    @debug "assignment made, starting greedy search"
    greedy_improve!(a; params = params)
    return a
end

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

function local_search!(a::Assignment, swap, params::GreedyParams)
    # select two nodes to swap and update data in the swap object
    make_swap!(swap, a, select_indices_swap(a, params.swap_rule))
    # apply swap, test if local improvement and update assignment if needed
    accept_reject_update!(a, swap, params.accept_rule)
end
