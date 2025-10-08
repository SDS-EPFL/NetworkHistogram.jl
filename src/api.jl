function nethist(data_input, dist_user, initial_node_labels,
        params::GreedyParams, zero_inflated::Bool = false)
    return _nethist(
        data_input, dist_user, initial_node_labels, params, Val(zero_inflated))
end

function _nethist(data_input, dist_user, initial_node_labels,
        params::GreedyParams, zero_inflated)
    @info "preprocessing data"
    dist = get_ref_dist(dist_user, zero_inflated)
    @show dist
    g = preprocess_data(data_input, dist, zero_inflated)

    @info "started optimization"
    out = greedy_optimize(g, initial_node_labels, params)

    @info "finished optimization with loglikelihood $(loglikelihood(out))"
    return postprocess(out)
end

function get_ref_dist(dist::D, ::Val{true}) where {D}
    return Dist(ZeroInflated(dist))
end
function get_ref_dist(dist::D, ::Val{false}) where {D}
    return Dist(dist)
end

function preprocess_data(data, dist::Dist, zero_inflated)
    A = EdgeList(_fast_compressed_obs(dist, data, zero_inflated))
    return A, dist
end

function postprocess(out)
    return out
    return out.node_labels, BlockModel(out)
end
