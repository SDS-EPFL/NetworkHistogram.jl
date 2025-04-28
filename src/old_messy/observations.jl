"""
    Observations{G, D}

A struct to hold observations for a network. The type parameter `G` represents the network
    structure and must support indexing and the `size` function.

# Fields
- `graph::G`: The network structure (e.g. adjacency matrix).
- `dist_ref::D`: distribution of the observations (used for getting support, type of elements, etc.)
"""
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

"""
    get_obs(graph::Observations, i::Int, j::Int)

Get the observation for the given pair of nodes.

# Arguments
- `graph::Observations`: The graph observations.
- `i::Int`: The first node.
- `j::Int`: The second node.

# Returns
- `obs`: The observation.
"""
function get_obs(graph::Observations, i::Int, j::Int)
    return graph.graph[i, j]
end

"""
    density(graph::Observations)

Get the density of the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `density`: The density of the graph.
"""
function density(graph::Observations)
    return sum(graph.graph) /
           ((size(graph.graph, 1) * (size(graph.graph, 1) - 1)))
end

"""
    get_degree(graph::Observations)

Get the degree of each node in the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `degrees`: The degrees of the nodes.
"""
function get_degree(graph::Observations)
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

"""
    normalized_laplacian(graph::Observations)

Get the normalized Laplacian of the graph.

# Arguments
- `graph::Observations`: The graph observations.

# Returns
- `L`: The normalized Laplacian matrix.
"""
function normalized_laplacian(graph::AbstractMatrix)
    degrees = sum(graph, dims = 1)
    degrees .-= minimum(degrees)
    n = size(graph, 1)
    L = similar(graph, Float64)
    for j in 1:n
        for i in 1:n
            if i == j
                L[i, j] = 1
            elseif degrees[i] == 0 || degrees[j] == 0
                L[i, j] = 0
            elseif graph[i, j] != 0
                L[i, j] = -1 / sqrt(degrees[i] * degrees[j])
            end
        end
    end
    return L
end

function Metis.graph(graph::Observations{
        G, <:UnivariateDistribution}) where {G}
    use_weights = true
    if minimum(graph.dist_ref) < 0
        @warn "Negative values are not allowed for MetisStart, using binary graph"
        use_weights = false
    end
    return Metis.graph(sparse(graph.graph), weights = use_weights)
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
function discretise(
        graph::Observations; number_groups = nothing, number_levels = nothing)
    if isnothing(number_groups) && isnothing(number_levels)
        throw(ArgumentError("Either `number_groups` or `number_levels` must be provided"))
    end
    if isnothing(number_levels)
        number_levels = round(Int,
            get_num_levels_from_groups(number_nodes(graph), number_groups))
    else
        if !isnothing(number_groups)
            @warn "disregarding `number_groups` as `number_levels` is provided"
        end
    end
    return discretise(
        graph, DiscretizerZeroToZero(number_levels, extrema(graph.graph)...))
end

"""
    discretise(graph::Observations, discretiser::Discretizer)

Discretise the graph observations using the given discretiser.

# Arguments
- `graph::Observations`: The graph observations.
- `discretiser::Discretizer`: The discretiser to use.

# Returns
- `discretised_graph`: The discretised graph observations.
- `discretiser`: The discretiser used.
"""
function discretise(graph::Observations, discretiser::Discretizer)
    A_encoded = encode(discretiser, graph.graph)
    return Observations(A_encoded, DiscretizedDistribution(discretiser)),
    discretiser
end

"""
    get_num_levels_from_groups(n, number_groups)

Get the number of levels for the discretized distribution given n and k.

# Arguments
- `n`: The number of nodes.
- `number_groups`: The number of groups.

# Returns
- `num_levels`: The number of levels.
"""
function get_num_levels_from_groups(n, number_groups)
    return max(1, n^(0.5 * (1 - log(number_groups) / log(n))))
end
