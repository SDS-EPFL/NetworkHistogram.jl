"""
    _default_init(dist::Distribution, start = MetisStart())

Initialize the distribution with a default rule.

# Arguments
- `dist::Distribution`: The distribution to initialize.
- `start`: The starting method.

# Returns
- `InitRule`: The initialization rule.
"""
function _default_init(dist::Distribution, start = MetisStart())
    if dist isa Bernoulli
        return InitRule(start, Val{BernoulliData}())
    elseif dist isa Categorical
        return InitRule(start, Val{CategoricalData}())
    elseif dist isa  DiscretizedDistribution || dist isa ZeroInflatedCategorical
        return InitRule(start, Val{SparseData}())
    else
        return InitRule(start, nothing)
    end
end

"""
    _nethist(g::Observations{G, D}, h; kwargs...)

Estimate the graphon and fit the model to the given graph observations.

# Arguments
- `g::Observations{G, D}`: The graph observations.
- `h`: Number of nodes per block.
- `kwargs...`: Additional keyword arguments.

# Returns
- `fit_model`: The fitted model.
- `a`: The assignment of nodes to blocks.
"""
function _nethist(g::Observations{G, D}, h; kwargs...) where {G, D}
    kwargs_dict = Dict(kwargs)
    start_clustering = pop!(kwargs_dict, :start_clustering, MetisStart())
    initialise_rule = pop!(
        kwargs_dict, :initialise_rule, _default_init(g.dist_ref, start_clustering))
    a = estimate_graphon(g, h;
        kwargs_dict..., initialise_rule = initialise_rule)
    return fit(a, g), a
end

"""
    nethist(g::Observations{G, D}; h, iterations, stalled_iter, swap_rule, accept_rule, progress_bar, start_clustering)

Fit a Stochastic Block Model (SBM) to the given graph observations.

# Arguments
- `g::Observations{G, D}`: The graph observations.
- `h`: Number of nodes per block.
- `iterations`: Maximum number of iterations.
- `stalled_iter`: Number of stalled iterations before stopping.
- `swap_rule::NodeSwapRule`: Rule for swapping nodes.
- `accept_rule::AcceptRule`: Rule for accepting swaps.
- `progress_bar::Bool`: Whether to show a progress bar.
- `start_clustering`: Initial clustering method.

# Returns
- `sbm`: The fitted SBM.
- `a`: The assignment of nodes to blocks.
"""
function nethist(g::Observations{G, D};
        h = select_number_node_per_block(g, EstimatedDegrees()),
        iterations = 100_000,
        stalled_iter = 1000,
        swap_rule::NodeSwapRule = RandomGroupSwap(),
        accept_rule::AcceptRule = Strict(),
        progress_bar::Bool = false,
        start_clustering = MetisStart()
) where {G, D}
    return _nethist(g, h;
        iterations = iterations,
        swap_rule = swap_rule,
        accept_rule = accept_rule,
        stop_rule = PreviousBestValue(stalled_iter),
        progress_bar = progress_bar,
        start_clustering = start_clustering)
end

"""
    nethist_discretised(g::Observations{G, D}; number_levels, h, iterations, stalled_iter, swap_rule, accept_rule, progress_bar, start_clustering)

Fit a discretised Stochastic Block Model (SBM) to the given graph observations.

# Arguments
- `g::Observations{G, D}`: The graph observations.
- `number_levels`: Number of levels for discretisation.
- `h`: Number of nodes per block.
- `iterations`: Maximum number of iterations.
- `stalled_iter`: Number of stalled iterations before stopping.
- `swap_rule::NodeSwapRule`: Rule for swapping nodes.
- `accept_rule::AcceptRule`: Rule for accepting swaps.
- `progress_bar::Bool`: Whether to show a progress bar.
- `start_clustering`: Initial clustering method.

# Returns
- `sbm_discretise`: The fitted discretised SBM.
- `a`: The assignment of nodes to blocks.
- `discretiser`: The discretiser used.
"""
function nethist_discretised(g::Observations{G, D};
        number_levels = nothing,
        h = select_number_node_per_block(g, EstimatedDegrees()),
        iterations = 100_000,
        stalled_iter = 1000,
        swap_rule::NodeSwapRule = RandomGroupSwap(),
        accept_rule::AcceptRule = Strict(),
        progress_bar::Bool = false,
        start_clustering = MetisStart()
) where {G, D}
    num_groups = isnothing(number_levels) ? number_nodes(g) ÷ h : nothing
    obs_discrete, discretiser = discretise(
        g, number_groups = num_groups, number_levels = number_levels)
    sbm_discretise, a = _nethist(obs_discrete, h;
        iterations = iterations,
        swap_rule = swap_rule,
        accept_rule = accept_rule,
        stop_rule = PreviousBestValue(stalled_iter),
        progress_bar = progress_bar,
        start_clustering = start_clustering)
    return sbm_discretise, a, discretiser
end
