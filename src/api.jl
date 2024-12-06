function _default_init(dist::Distribution, start = MetisStart())
    if dist isa Bernoulli
        return InitRule(start, Val{BernoulliData}())
    elseif dist isa Categorical || dist isa CategoricalArray ||
           dist isa DiscretizedDistribution
        return InitRule(start, Val{CategoricalData}())
    else
        return InitRule(start, nothing)
    end
end

function _nethist(g::Observations{G, D}, h; kwargs...) where {G, D}
    kwargs_dict = Dict(kwargs)
    start_clustering = pop!(kwargs_dict, :start_clustering, MetisStart())
    initialise_rule = pop!(
        kwargs_dict, :initialise_rule, _default_init(g.dist_ref, start_clustering))
    a = estimate_graphon(g, h;
        kwargs_dict..., initialise_rule = initialise_rule)
    return fit(a, g), a
end

function nethist(g::Observations{G, D};
        h = select_number_node_per_block(g, EstimatedDegrees()),
        max_iter = 100_000,
        stalled_iter = 1000,
        swap_rule::NodeSwapRule = RandomGroupSwap(),
        accept_rule::AcceptRule = Strict(),
        progress_bar::Bool = false,
        start_clustering = MetisStart()
) where {G, D}
    return _nethist(g, h;
        max_iter = max_iter,
        swap_rule = swap_rule,
        accept_rule = accept_rule,
        stop_rule = PreviousBestValue(stalled_iter),
        progress_bar = progress_bar,
        start_clustering = start_clustering)
end

function nethist_discretised(g::Observations{G, D};
        number_levels = nothing,
        h = select_number_node_per_block(g, EstimatedDegrees()),
        max_iter = 100_000,
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
        max_iter = max_iter,
        swap_rule = swap_rule,
        accept_rule = accept_rule,
        stop_rule = PreviousBestValue(stalled_iter),
        progress_bar = progress_bar,
        start_clustering = start_clustering)
    return sbm_discretise, a, discretiser
end
