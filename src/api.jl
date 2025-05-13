function nethist(data_input, dist_user, initial_node_labels, params::GreedyParams)

    dist = Dist(dist_user)
    g = preprocess_data(data_input, dist)

    @debug "started optimizatiion"
    out = greedy_optimize(g, initial_node_labels, params)

    return postprocess(out)
end


function preprocess_data(data, dist::Dist)
    A = EdgeList(_fast_compressed_obs(dist, data))
    return  A, dist
end


function postprocess(out)
    return out
    return BlockModel(out)
end
