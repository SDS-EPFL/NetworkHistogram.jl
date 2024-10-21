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

function normalized_laplacian(g::Observations)
    return normalized_laplacian(g.graph)
end

function normalized_laplacian(g::AbstractGraph)
    return normalized_laplacian(Graphs.adjacency_matrix(g))
end

function normalized_laplacian(g::AbstractMatrix)
    degrees = sum(g, dims = 1)
    n = size(g, 1)
    L = similar(g, Float64)
    for j in 1:n
        for i in 1:n
            if i == j
                L[i, j] = 1
            elseif g[i, j] == 1
                L[i, j] = -1 / sqrt(degrees[i] * degrees[j])
            end
        end
    end
    return L
end
