function nethist(data_input, dist_user, initial_node_labels, params::GreedyParams, zero_inflated::Bool = false)
    return _nethist(data_input, dist_user, initial_node_labels, params, Val(zero_inflated))
end


function _nethist(data_input, dist_user, initial_node_labels, params::GreedyParams, zero_inflated)
    @debug "preprocessing data"
    dist = get_ref_dist(dist_user, zero_inflated)
    g = preprocess_data(data_input, dist)

    @debug "started optimizatiion"
    out = greedy_optimize(g, initial_node_labels, params)

    @debug "finished optimizatiion with loglikelihood $(loglikelihood(out))"
    return postprocess(out)
end


function get_ref_dist(dist::D, ::Val{true}) where {D}
    return Dist(ZeroInflated(dist))
end
function get_ref_dist(dist::D, ::Val{false}) where {D}
    return Dist(dist)
end

function preprocess_data(data, dist::Dist)
    A = EdgeList(_fast_compressed_obs(dist, data))
    return  A, dist
end


function postprocess(out)
    return out
    return out.node_labels, BlockModel(out)
end
