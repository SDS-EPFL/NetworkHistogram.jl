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

using NetworkHistogram, Random

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
obs = NetworkHistogram.Observations(G, dist)
sbm_fit = NetworkHistogram.fit(a, obs)

sbm = NetworkHistogram.initialize_sbm([1 / 3, 1 / 3, 1 / 3], dist)
for i in 1:3
    sbm[i, i] = Bernoulli(0.2)
    for j in (i + 1):3
        sbm[i, j] = Bernoulli(0.01)
    end
end

size_per_block = 200
A, node_labels = NetworkHistogram.sample(sbm, 3 * size_per_block);
node_labels = repeat(1:3, inner = size_per_block)
Random.shuffle!(node_labels)
group_number = NetworkHistogram.GroupSize(size(A, 1), size_per_block)
a_star = NetworkHistogram.Assignment(group_number, node_labels, additional_info)
obs_star = NetworkHistogram.Observations(SimpleGraph(A), dist)
sbm_fitted = NetworkHistogram.fit(a_star, obs_star)

sbm_fitted

init_rule = NetworkHistogram.InitRule(NetworkHistogram.RandomStart(), nothing)
ll_old = NetworkHistogram.score(a_star, obs_star)
println("Log likelihood: ", ll_old)
a_best = NetworkHistogram.optimize(obs_star, size_per_block, max_iter = 100;
    stop_rule = NetworkHistogram.PreviousBestValue(5),
    initialise_rule = init_rule, progress_bar = true)
ll_new = NetworkHistogram.score(a_best, obs_star)
println("Log likelihood: ", ll_new)
if ll_old < ll_new
    println("Optimization improved the log likelihood.")
else
    println("Optimization did not improve the log likelihood.")
end
