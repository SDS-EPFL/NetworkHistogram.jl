#=
# Decorated Graphon Tutorial for Weighted Networks
=#
using Clustering
using NetworkHistogram
using Distributions
using LinearAlgebra
using Random

import Distributions: pdf

pdf_kuma(α, β, x, p = 1.0) = @. p * (α * β * x^(α - 1) .* (1 - x^α)^(β - 1))

graphon_params = (x, y) -> (4 * (cos(π * (x - y)) + 1) + 1, max(x, y) * 8 + 1)

graphon = DecoratedGraphon((x, y) -> Kumaraswamy(graphon_params(x, y)...))

import CairoMakie as Mke
let
    fig = Mke.Figure()
    ax = Mke.Axis(fig[1, 1], aspect = Mke.DataAspect())
    hm = Mke.heatmap!(ax, graphon, colormap = :viridis)
    Mke.Colorbar(fig[1, 2], hm)
    ax2 = Mke.Axis(fig[1, 3], aspect = Mke.DataAspect())
    hm2 = Mke.heatmap!(ax2, graphon, k = 2, colormap = :viridis)
    Mke.Colorbar(fig[1, 4], hm2)
    fig
end

n = 2000
k = 5
n_bins = 20
p = 0.9
A = sample_graph(graphon, n) .* Symmetric(rand(Bernoulli(p), n, n));
ξs = range(0, 1; length = n)
oracle_latents = ordered_start_labels(n, k);

res_oracle = NetworkHistogram.oracle_estimator(
    A, oracle_latents, NetworkHistogram.UnitIntervalConvertor(n_bins));

starting_labels = shuffle(oracle_latents);

max_iter = 1_000_000
stalled_iters = 5_000

res_new = NetworkHistogram.nethist_continuous(
    A, k,
    starting_labels;
    bins = n_bins
);

# convertor = NetworkHistogram.UnitIntervalConvertor(10)

# data = convertor.(A)
# es_new = NetworkHistogram.GreedySuffStats(
#     data, initial_labels, num_categories = num_bins(convertor),
#     type_suff_stats = :categorical,
#     max_iter = max_iter,
#     swap_rule = NetworkHistogram.RandomGroupSwap(),
#     stop_rule = NetworkHistogram.PreviousBestValue(stalled_iters, Inf, :min),
#     progress = true
# );
# node_labels_es_new, parameters = NetworkHistogram.estimate!(
#     es_new, data, initial_labels; iter_progress = 10_000)

# model_es_new = NetworkHistogram.DecoratedSBM(to_distribution.(convertor, parameters),
#     counts(node_labels_es_new) ./ length(node_labels_es_new));

# res_new = NetworkHistogram.NethistResult(node_labels_es_new, model_es_new);
NetworkHistogram.align_res_true_latents!(res_new, res_oracle.labels);
xs = range(0, 1; length = 100)

function viz_one_group!(axis, g1, g2, A, ξs, res_oracle, res_new, xs; n_viz = 20, p = p)
    nodes_1 = findall(res_oracle.labels .== g1)
    nodes_2 = findall(res_oracle.labels .== g2)
    edge_values = [A[x, y] for y in nodes_2 for x in nodes_1]
    Mke.vlines!(axis, edge_values, ymax = 0.025, color = :lightgray)
    # Mke.hist!(axis, edge_values; normalization = :pdf, color = :gray)
    x1 = sample(ξs[nodes_1], n_viz, replace = false)
    x2 = sample(ξs[nodes_2], n_viz, replace = false)
    for x_ in x1
        for y_ in x2
            Mke.lines!(axis, xs, pdf_kuma(graphon_params(x_, y_)..., xs, p),
                color = :gray, alpha = 0.1)
        end
    end
    Mke.lines!(axis, xs, map(Base.Fix1(pdf, res_oracle.model.θ[g1, g2]), xs),
        color = :blue, label = "True")
    Mke.lines!(axis, xs, map(Base.Fix1(pdf, res_new.model.θ[g1, g2]), xs),
        color = :black, linestyle = :dash, label = "Estimated")
end

for g in 1:k
    for g2 in 1:g
        fig = Mke.Figure(size = (600, 400))
        ax = Mke.Axis(fig[1, 1], title = "Group $g vs Group $g2", xlabel = "Edge Value",
            ylabel = "Density")
        viz_one_group!(ax, g, g2, A, ξs, res_oracle, res_new, xs, p = p, n_viz = 5)
        display(fig)
    end
end

##
ssm_test = SSM(res_new.model, k)

shape_range = 1:(k * (k + 1) ÷ 2 - 1)
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res_new.model, A, res_new.labels, shape_range)

Mke.lines(shape_range, criterion_values)

##
# using Kneedle
# kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing)
# #  Let's extract the optimal number of shapes using the Kneedle algorithm:

# k_knee = knees(kr)[1]
# ssm_knee = SSM(res_new.model, k_knee)

##

clustering_res = kmeans(A, k)

res_kmeans = NetworkHistogram.oracle_estimator(
    A, assignments(clustering_res), NetworkHistogram.UnitIntervalConvertor(n_bins);
    type_suff_stats = Val(:categorical),
    name = "k-means");

NetworkHistogram.align_res_true_latents!(res_kmeans, res_oracle.labels);

for g in 1:k
    for g2 in 1:g
        fig = Mke.Figure(size = (600, 400))
        ax = Mke.Axis(fig[1, 1], title = "Group $g vs Group $g2", xlabel = "Edge Value",
            ylabel = "Density")
        viz_one_group!(ax, g, g2, A, ξs, res_oracle, res_kmeans, xs, p = p, n_viz = 5)
        display(fig)
    end
end
