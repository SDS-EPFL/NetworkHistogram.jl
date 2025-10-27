#=
# Decorated Graphon Tutorial for Weighted Networks
=#
using Clustering
using NetworkHistogram
using Distributions
using LinearAlgebra
using Random

graphon = DecoratedGraphon((x, y) -> Kumaraswamy(
    4 * (cos(π * (x - y)) + 1) + 1, max(x, y) * 8 + 1))

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

n = 4000
k = 4
A = sample_graph(graphon, n) .* Symmetric(rand(Bernoulli(0.9), n, n));
oracle_latents = ordered_start_labels(n, k);
starting_labels = copy(oracle_latents);
p_shuffle = 1 - 1.5 / k
@info "Shuffling $(p_shuffle*100)% labels for starting point"
indices_to_shuffle = sample(1:n, floor(Int, n * p_shuffle), replace = false);
starting_labels[indices_to_shuffle] .= shuffle(starting_labels[indices_to_shuffle]);
@assert starting_labels != oracle_latents

max_iter = 1_000_000
stalled_iters = 5_000

res_new = NetworkHistogram.nethist_continuous(
    A, k,
    starting_labels;
    num_bins_ = 10,
    max_iter = max_iter,
    stalled_iters = stalled_iters,
    progress_bar = true
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
NetworkHistogram.align_res_true_latents!(res_new, oracle_latents);

##
ssm_test = SSM(res_new.model, k)

shape_range = 1:min(5, k * (k + 1) ÷ 2 - 1)
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res_new.model, A, res_new.labels, shape_range)

Mke.lines(shape_range, criterion_values)

##
# using Kneedle
# kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing)
# #  Let's extract the optimal number of shapes using the Kneedle algorithm:

# k_knee = knees(kr)[1]
# ssm_knee = SSM(res_new.model, k_knee)
