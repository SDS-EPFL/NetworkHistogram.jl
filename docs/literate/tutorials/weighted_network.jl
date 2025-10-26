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

n = 5000
k = 15
A = sample_graph(graphon, n) .* Symmetric(rand(Bernoulli(0.9), n, n));
oracle_latents = ordered_start_labels(n, k);
starting_labels = copy(oracle_latents);
# shuffle!(starting_labels);
p_shuffle = 1 - 1.5 / k
@info "Shuffling $(p_shuffle*100)% labels for starting point"
indices_to_shuffle = sample(1:n, floor(Int, n * p_shuffle), replace = false);
starting_labels[indices_to_shuffle] .= shuffle(starting_labels[indices_to_shuffle]);
@assert starting_labels != oracle_latents

res, res_cat, A_cat = NetworkHistogram.nethist_continuous_edges(A,
    starting_labels, GreedyParams(
        1_000_000,
        RandomGroupSwap(),
        Strict(),
        PreviousBestValue(10_000, Inf, :min),
        true # progress bar
    );
    num_bins_ = 10, lower_bound = eps(), upper_bound = 1);

latents = range(0, 1; length = n);

ssm_test = SSM(res.model, k)

shape_range = 1:min(30, k * (k + 1) ÷ 2 - 1)
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res_cat.model, A_cat, latents, shape_range)

Mke.lines(shape_range, criterion_values)

##
using Kneedle
kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing)
#  Let's extract the optimal number of shapes using the Kneedle algorithm:

k_knee = knees(kr)[1]
ssm_knee = SSM(res.model, k_knee)
