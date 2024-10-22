# Slow fallback methods for the Assignment type
# speed up by implementing specialized methods for the BernoulliAssignment type and others
# method to compute estimator from node clustering as specified in assignment
function fit(a::Assignment, g::Observations)
    dists = initialize_sbm(a.group_size, g.dist_ref)
    for group1 in 1:number_groups(a)
        for group2 in group1:number_groups(a)
            edge_indices = get_edge_indices(a, group1, group2)
            dists[group1,
            group2] = fit_group(g.dist_ref, g, edge_indices)
        end
    end
    return dists
end

function fit_group(distribution, g, edges)
    return Distributions.fit(typeof(distribution), get_obs.(Ref(g), edges))
end

# method to compute the log likelihood of a SBM fitted according to the assignment
function loglikelihood(a::Assignment, g::Observations)
    return _log_likelihood(a, fit(a, g), g)
end

function _log_likelihood(a::Assignment, sbm::SBM, g)
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
