#
struct Observations{G, D}
    graph::G
    dist_ref::D
end

function number_nodes(g::Observations{AbstractGraph, D}) where {D}
    return nv(g.graph)
end

function number_nodes(g::Observations)
    return size(g.graph,1)
end
function get_obs(g::Observations, x::Tuple)
    return get_obs(g.graph, x[1], x[2])
end

function get_obs(g::SimpleGraph, x::Tuple)
    return get_obs(g, x[1], x[2])
end

function get_obs(g::SimpleGraph, i::Int, j::Int)
    return convert(Bool, has_edge(g, i, j))
end

get_obs(g::AbstractArray, x) = get_obs(g, x[1], x[2])
get_obs(g::AbstractArray, i, j) = g[i, j]
