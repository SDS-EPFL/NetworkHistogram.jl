# remove all references to graphs, and only use sparse matrices ?
struct Observations{G, D}
    graph::G
    dist_ref::D
end

"""
    number_nodes(graph::Observations)

Get the number of nodes in the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `num_nodes`: The number of nodes.
"""
function number_nodes(graph::Observations{AbstractGraph, D}) where {D}
    return nv(graph.graph)
end

function number_nodes(graph::Observations)
    return size(graph.graph, 1)
end

"""
    get_obs(graph::Observations, x::Tuple)

Get the observation for the given tuple of nodes.

# Arguments
- `graph::Observations`: The graph observations.
- `x::Tuple`: The tuple of nodes.

# Returns
- `obs`: The observation.
"""
function get_obs(graph::Observations, x::Tuple)
    return get_obs(graph, x[1], x[2])
end

function get_obs(graph::Observations, i::Int, j::Int)
    return get_obs(graph.graph, i, j)
end

function get_obs(g::SimpleGraph, x::Tuple)
    return get_obs(g, x[1], x[2])
end

function get_obs(g::SimpleGraph, i::Int, j::Int)
    return convert(Bool, has_edge(g, i, j))
end

get_obs(g::AbstractArray, x) = get_obs(g, x[1], x[2])
get_obs(g::AbstractArray, i, j) = g[i, j]

"""
    density(graph::Observations)

Get the density of the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `density`: The density of the graph.
"""
density(graph::Observations) = density(graph.graph)
function density(g::AbstractGraph)
    return Graphs.density(g)
end

function density(g::AbstractMatrix)
    return sum(g) / ((size(g, 1) * (size(g, 1) - 1)))
end

"""
    get_degree(graph::Observations)

Get the degree of each node in the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `degrees`: The degrees of the nodes.
"""
function get_degree(graph::Observations{AbstractGraph, D}) where {D}
    Graphs.degree(graph.graph)
end

function get_degree(graph)
    return sum(graph.graph, dims = 2)
end

"""
    get_adj(graph::Observations)

Get the adjacency matrix of the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `adj_matrix`: The adjacency matrix.
"""
function get_adj(graph::Observations{AbstractGraph, D}) where {D}
    return Graphs.adjacency_matrix(graph.graph)
end

function get_adj(graph::Observations)
    return graph.graph
end

function normalized_laplacian(graph::Observations)
    return normalized_laplacian(graph.graph)
end

function normalized_laplacian(g::AbstractGraph)
    return normalized_laplacian(Graphs.adjacency_matrix(g))
end


normalized_laplacian(g::CategoricalArray) = normalized_laplacian(levelcode.(g))

function normalized_laplacian(g::AbstractMatrix)
    degrees = sum(g, dims = 1)
    degrees .-= minimum(degrees)
    n = size(g, 1)
    L = similar(g, Float64)
    for j in 1:n
        for i in 1:n
            if i == j
                L[i, j] = 1
            elseif degrees[i] == 0 || degrees[j] == 0
                L[i, j] = 0
            elseif g[i, j] != 0
                L[i, j] = -1 / sqrt(degrees[i] * degrees[j])
            end
        end
    end
    return L
end

function Metis.graph(graph::Observations{<:AbstractGraph, <:Bernoulli})
    return Metis.graph(graph.graph)
end

function Metis.graph(g::Observations{<:AbstractMatrix, <:Bernoulli})
    return Metis.graph(SimpleGraph(g.graph))
end

function Metis.graph(graph::Observations{<:AbstractGraph, <:UnivariateDistribution})
    if minimum(graph.dist_ref) < 0
        @warn "Negative values are not allowed for MetisStart, using binary graph"
        return Metis.graph(graph.graph)
    else
        return Metis.graph(graph.graph, weights = true)
    end
end

function Metis.graph(g::Observations{<:AbstractMatrix, <:UnivariateDistribution})
    return Metis.graph(
        weights(SimpleWeightedGraph(g.graph)), weights = true)
end

function Metis.graph(g::Observations{<:CategoricalMatrix, <:UnivariateFinite})
    A, _ = categorical_matrix(g)
    return Metis.graph(
        adjacency_matrix(SimpleWeightedGraph(A)), weights = true)
end


"""
    discretise(graph::Observations; number_groups, number_levels)

Discretise the graph observations.

# Arguments
- `graph::Observations`: The graph observations.
- `number_groups`: Number of groups for discretisation.
- `number_levels`: Number of levels for discretisation.

# Returns
- `discretised_graph`: The discretised graph observations.
- `discretiser`: The discretiser used.


Assume that the diagonal is zero.
0 indicates no edge, while missing indicates no information about the edge.
By default maps 0 to 0. If you want another behaviour use the function where you
pass a `Discretizer` object.

number_levels will be the number of levels in the discretized distribution (excluding 0).
"""
function discretise(graph::Observations{G, D};
        number_groups = nothing, number_levels = nothing) where {G, D}
    if isnothing(number_groups) && isnothing(number_levels)
        throw(ArgumentError("Either `number_groups` or `number_levels` must be provided"))
    end
    if isnothing(number_levels)
        number_levels = round(Int,get_num_levels_from_groups(number_nodes(graph), number_groups))
    else
        if !isnothing(number_groups)
            @warn "disregarding `number_groups` as `number_levels` is provided"
        end
    end
    return discretise(graph, DiscretizerZeroToZero(number_levels, extrema(graph.graph)...))
end

function discretise(graph::Observations{G, D}, discretiser ::Discretizer) where {G,D<:UnivariateDistribution}
    A_encoded = encode(discretiser, _graph_to_mat(graph))
    return Observations(A_encoded, DiscretizedDistribution(discretiser)), discretiser
end


function _graph_to_mat(graph::Observations{<:AbstractGraph, D}) where {D<:UnivariateDistribution}
    return weights(graph.graph)
end

function _graph_to_mat(graph::Observations{<:AbstractMatrix, D}) where {D<:UnivariateDistribution}
    return graph.graph
end


"""
Get the number of levels for the discretized distribution given n and k.
"""
function get_num_levels_from_groups(n, number_groups)
    return max(1,  n^(0.5 * (1 - log(number_groups) / log(n))))
end
