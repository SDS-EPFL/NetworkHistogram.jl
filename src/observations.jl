#
struct Observations{G, D}
    graph::G
    dist_ref::D
end

function number_nodes(g::Observations{AbstractGraph, D}) where {D}
    return nv(g.graph)
end

function number_nodes(g::Observations)
    return size(g.graph, 1)
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

density(g::Observations) = density(g.graph)
function density(g::AbstractGraph)
    return Graphs.density(g)
end

function density(g::AbstractMatrix{Bool})
    return sum(g) / ((size(g, 1) * (size(g, 1) - 1)))
end

function get_degree(g::Observations{AbstractGraph, D}) where {D}
    Graphs.degree(g.graph)
end

function get_degree(g)
    return sum(g.graph, dims = 2)
end

function get_adj(g::Observations{AbstractGraph, D}) where {D}
    return Graphs.adjacency_matrix(g.graph)
end

function get_adj(g::Observations)
    return g.graph
end
