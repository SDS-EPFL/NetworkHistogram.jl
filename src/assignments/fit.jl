# method to compute estimator from node clustering as specified in assignment
function fit(a::Assignment, g::Observations)
    dists = initialize_sbm(a.group_size, g.dist_ref)
    for group1 in 1:number_groups(a)
        for group2 in group1:number_groups(a)
            edge_indices = get_edge_indices(a, group1, group2)
            dists[group1, group2] = fit(g.dist_ref, g.graph, edge_indices)
        end
    end
    return dists
end

function fit(distribution, g, edges)
    return Distributions.fit(typeof(distribution), get_obs.(Ref(g), edges))
end

# method to compute the log likelihood of fitted SBM
function log_likelihood(a::Assignment, sbm::SBM, g)
    log_likelihood = 0.0
    for i in 1:number_nodes(a)
        label_a = a.node_labels[i]
        for j in (i + 1):number_nodes(a)
            label_b = a.node_labels[j]
            log_likelihood += logdensityof(sbm[label_a, label_b], get_obs(g, i, j))
        end
    end
    return log_likelihood
end

function log_likelihood(a::Assignment, g::Observations)
    dists = fit(a, g)
    return log_likelihood(a, dists, g.graph)
end

# default score is the log likelihood
function score(a::Assignment, g::Observations)
    return log_likelihood(a, g)/binomial(number_nodes(a), 2)
end
