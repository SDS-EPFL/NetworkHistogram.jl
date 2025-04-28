function nethist(data_input, dist_user, initial_node_labels, params::GreedyParams)

    dist = Dist(dist_user)
    g = preprocess_data(data_input, dist)


    out = greedy_optimize(g, initial_node_labels, params)

    return postprocess(out)
end


function preprocess_data(data, dist)
    A = _fast_compressed_g.(dist, data)
    return  A, dist
end


function postprocess(out)
    return true
    return BlockModel(optimal_a)
end
