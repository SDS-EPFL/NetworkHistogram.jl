using Graphs
using MetaGraphsNext
using Distributions
using SimpleWeightedGraphs

cities = MetaGraph(
    Graph();
    label_type = String,
    vertex_data_type = Vector{Float64},
    edge_data_type = Float64,
    graph_data = nothing,
    default_weight = -Inf,
    weight_function = identity,
);

cities["Paris"] = [2.3, 48.9];
cities["London"] = [0.1, 51.5];
cities["Berlin"] = [13.4, 52.5];
cities["Lausanne"] = [6.6, 46.5];
cities["Paris", "London"] = 0.5;
cities["Paris", "Berlin"] = 1;
cities["London", "Berlin"] = 0;
cities["Paris", "Lausanne"] = 1;

getindex(cities, "Paris")


typeof(cities)

index = code_for(cities, "Paris")
label_for.(Ref(cities), neighbors(cities, index))
adjacency_matrix(cities)
collect(weights(cities))

G = Graph(20, 20)
adjacency_matrix(G)


sources = [1, 2, 1];

destinations = [2, 3, 3];

weight_edges = [0.5, 0.8, 2.0];

g = SimpleWeightedGraph(sources, destinations, weight_edges)
add_vertices!(g, 5)


function get_obs(g::SimpleGraph{T}, i::Int, j::Int) where {T}
    return has_edge(g, i, j)
end
