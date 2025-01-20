# Slow fallback methods for the Assignment type
# Speed up by implementing specialized methods for the BernoulliAssignment type and others

"""
    fit(a::Assignment, g::Observations)

Compute the estimator from node clustering as specified in the assignment.

# Arguments
- `a::Assignment`: The assignment of nodes to blocks.
- `g::Observations`: The graph observations.

# Returns
- `dists`: The fitted distributions.
"""
function fit(a::Assignment, g::Observations)
    dists = initialize_sbm(a.group_size, g.dist_ref)
    fit!(dists, g, a)
    return dists
end

"""
    fit!(sbm::BlockModel{D,K,F}, g::Observations{G,D}, a::Assignment) where {G,D,K,F}

Fit the SBM to the given graph observations and assignment.

# Arguments
- `sbm::BlockModel{D,K,F}`: The block model to fit.
- `g::Observations{G,D}`: The graph observations.
- `a::Assignment`: The assignment of nodes to blocks.
"""
function fit!(sbm::BlockModel{D, K, F}, g::Observations{G, D},
        a::Assignment) where {G, D, K, F}
    for group1 in 1:number_groups(a)
        for group2 in group1:number_groups(a)
            edge_indices = get_edge_indices(a, group1, group2)
            sbm[group1, group2] = fit_group(g.dist_ref, g, edge_indices)
        end
    end
end

function fit_group(d::ZeroInflatedCategorical, g, edges)
    return Distributions.fit(
        typeof(d), get_obs.(Ref(g), edges), ncategories(g.dist_ref))
end

function fit_group(distribution, g, edges)
    return Distributions.fit(typeof(distribution), get_obs.(Ref(g), edges))
end

function fit_group(distribution::Binomial, g, edges)
    return Distributions.fit(
        typeof(distribution), ntrials(distribution), get_obs.(Ref(g), edges))
end

"""
    loglikelihood(a::Assignment, g::Observations)

Compute the log likelihood of a BlockModel fitted according to the assignment.

# Arguments
- `a::Assignment`: The assignment of nodes to blocks.
- `g::Observations`: The graph observations.

# Returns
- `log_likelihood`: The log likelihood of the fitted model.
"""
function loglikelihood(a::Assignment, g::Observations)
    return _log_likelihood(a, fit(a, g), g)
end

function _log_likelihood(a::Assignment, sbm::BlockModel, g)
    log_likelihood = 0.0
    for i in 1:number_nodes(a)
        label_a = a.node_labels[i]
        for j in (i + 1):number_nodes(a)
            label_b = a.node_labels[j]
            log_likelihood += logdensityof(
                sbm[label_a, label_b], get_obs(g, i, j))
        end
    end
    return log_likelihood
end

"""
    fit!(sbm::BlockModel{D,K,F}, g::Observations{G,D}) where {G,D,K,F}

Fit the SBM to the given graph observations.

# Arguments
- `sbm::BlockModel{D,K,F}`: The block model to fit.
- `g::Observations{G,D}`: The graph observations.
"""
function fit!(
        sbm::BlockModel{D, K, F}, g::Observations{G, D}) where {G, D, K, F}
    k = number_blocks(sbm)
    a = estimate_graphon(g, select_number_node_per_block(g, OracleK(k)))
    fit!(sbm, g, a)
end
