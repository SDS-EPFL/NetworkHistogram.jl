#=
# Decorated Graphon Tutorial for Multiplex Networks
=#
using NetworkHistogram
using Distributions
using StaticArrays
import CairoMakie as Mke

using Random
Random.seed!(1234);
h = 300;

function W_multiplex(x, y)
    ps = zeros(4)
    ps[2] = sqrt(abs(x - y)) / 2           # layer 1 only
    ps[3] = abs(sin(2π * x) * sin(2π * y)) / 2  # layer 2 only
    ps[4] = min(x, y) / 4                   # both layers
    ps[1] = 1 - sum(ps[2:4])                # no edge
    return DiscreteNonParametric(0:3, SVector{4}(ps))
end

graphon = DecoratedGraphon(W_multiplex)

let
    fig = Mke.Figure(size = (4 * h, h))
    for m in 1:4
        ax = Mke.Axis(fig[1, m], aspect = Mke.DataAspect())
        Mke.heatmap!(ax, graphon, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
end

n = 1000
true_latents = range(0, 1; length = n)
A = sample_graph(graphon, true_latents);

k = 20
oracle_labels = ordered_start_labels(n, k);
initial_labels = shuffle(oracle_labels);

res = NetworkHistogram.nethist_discrete_edges(A,
    initial_labels, GreedyParams(
        1_000_000,
        RandomGroupSwap(),
        Strict(),
        PreviousBestValue(5_000, Inf, :min),
        true
    ));

let
    fig = Mke.Figure(size = (4 * h, h))
    for m in 1:4
        ax = Mke.Axis(fig[1, m], aspect = Mke.DataAspect())
        Mke.heatmap!(ax, res.model, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
end

NetworkHistogram.align_res_true_latents!(res, oracle_labels);
let
    fig = Mke.Figure(size = (4 * h, h))
    for m in 1:4
        ax = Mke.Axis(fig[1, m], aspect = Mke.DataAspect())
        Mke.heatmap!(ax, res.model, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
end

using Clustering
shape_range = 1:20
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res.model, A, true_latents, shape_range);

using Kneedle
# kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing)
# #  Let's extract the optimal number of shapes using the Kneedle algorithm:

# k_knee = knees(kr)[1]
# k_knee = 10
ssm = SSM(res.model, k_knee)

let
    fig = Mke.Figure(size = (4 * h, 3 * h))
    for (i, model) in enumerate([graphon, res.model, ssm])
        for m in 1:4
            ax = Mke.Axis(fig[i, m], aspect = Mke.DataAspect())
            Mke.heatmap!(ax, model, k = m, colormap = :binary, colorrange = (0, 1))
        end
    end
    fig
end
