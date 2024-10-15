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
    weight_function = identity
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

##

using NetworkHistogram
group_number = NetworkHistogram.GroupSize(nv(G), 3)
if typeof(group_number) == NetworkHistogram.GroupSize{Tuple{Int, Int}}
    node_labels = repeat(1:(length(group_number) - 1), inner = group_number[1])
    last_labels = fill(length(group_number), group_number[end])
    node_labels = vcat(node_labels, last_labels)
else
    node_labels = repeat(1:length(group_number), inner = group_number[1])
end
additional_info = 1
a = NetworkHistogram.Assignment(group_number, node_labels)
dist = Bernoulli(0.5)
sbm_fit = NetworkHistogram.fit(a, G, dist)

sbm = NetworkHistogram.initialize_sbm([1 / 3, 1 / 3, 1 / 3], dist)
for i in 1:3
    sbm[i, i] = Bernoulli(0.8)
    for j in (i + 1):3
        sbm[i, j] = Bernoulli(0.01 + 0.1 * (i + j))
    end
end

size_per_block = 200
A, node_labels = NetworkHistogram.sample(sbm, 3 * size_per_block);
node_labels = repeat(1:3, inner = size_per_block)
group_number = NetworkHistogram.GroupSize(size(A, 1), size_per_block)
a_star = NetworkHistogram.Assignment(group_number, node_labels, additional_info)
sbm_fitted = NetworkHistogram.fit(a_star, SimpleGraph(A), dist)

sbm_fitted
