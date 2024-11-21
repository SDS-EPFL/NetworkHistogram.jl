# switch to MetaGraphsNext.jl ?
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
    return get_obs(g, x[1], x[2])
end

function get_obs(g::Observations, i::Int, j::Int)
    return get_obs(g.graph, i, j)
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

function density(g::AbstractMatrix)
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

function Metis.graph(g::Observations{<:AbstractGraph, <:Bernoulli})
    return Metis.graph(g.graph)
end

function Metis.graph(g::Observations{<:AbstractMatrix, <:Bernoulli})
    return Metis.graph(SimpleGraph(g.graph))
end

function Metis.graph(g::Observations{<:AbstractMatrix, <:Categorical})
    return Metis.graph(
        adjacency_matrix(SimpleWeightedGraph(g.graph)), weights = true)
end

function Metis.graph(g::Observations{<:CategoricalMatrix, <:UnivariateFinite})
    A, _ = categorical_matrix(g)
    return Metis.graph(
        adjacency_matrix(SimpleWeightedGraph(A)), weights = true)
end

function discretise(g::Observations{<:AbstractMatrix{R}, D};
        number_groups = nothing, number_levels = nothing) where {R<:Real, D}
    if isnothing(number_groups) && isnothing(number_levels)
        throw(ArgumentError("Either `number_groups` or `number_levels` must be provided"))
    end
    if isnothing(number_levels)
        number_levels = round(Int,get_num_levels_from_groups(number_nodes(g), number_groups))
    else
        if !isnothing(number_groups)
            @warn "disregarding `number_groups` as `number_levels` is provided"
        end
    end
    #zero_locations = g.graph .== 0
    bin_edges = binedges(DiscretizeUniformWidth(number_levels), g.graph)
    discretizer = LinearDiscretizer(bin_edges)
    A_encoded = encode(discretizer, g.graph)
    for i in 1:size(A_encoded, 1)
        A_encoded[i, i] = 0
    end
    #A_encoded[zero_locations] .= 0
    return Observations(A_encoded, Categorical(number_levels + 1)), discretizer
end

function get_num_levels_from_groups(n, number_groups)
    return n^(0.5 * (1 - log(number_groups) / log(n)))
end
