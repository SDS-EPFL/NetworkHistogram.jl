#=
# Decorated Graphon Tutorial for Weighted Networks
=#
using Clustering
using NetworkHistogram
using Distributions
using LinearAlgebra
using Random
using Graphons

import Distributions: pdf

pdf_kuma(α, β, x, p = 1.0) = @. p * (α * β * x^(α - 1) .* (1 - x^α)^(β - 1))

graphon_params = (x, y) -> (3 * abs(sin(2 * π * x) * sin(2 * π * y)) + 0.8, max(x, y) * 8)

graphon = DecoratedGraphon((x, y) -> Kumaraswamy(graphon_params(x, y)...))

import CairoMakie as Mke
let
    fig = Mke.Figure(size = (510, 200))
    ax = Mke.Axis(fig[1, 1], aspect = Mke.DataAspect(), title = "α")
    hm = Mke.heatmap!(ax, graphon, k = 1, colormap = :viridis)
    Mke.Colorbar(fig[1, 2], hm)
    ax2 = Mke.Axis(fig[1, 3], aspect = Mke.DataAspect(), title = "β")
    hm2 = Mke.heatmap!(ax2, graphon, k = 2, colormap = :viridis)
    Mke.Colorbar(fig[1, 4], hm2)
    fig
end

# We sample a weighted network from the graphon

Random.seed!(1234);
n = 2000
k = 12
n_bins = 20
p = 0.8

A = sample_graph(graphon, n) .* Symmetric(rand(Bernoulli(p), n, n));
ξs = range(0, 1; length = n)
oracle_latents = ordered_start_labels(n, k);

res_oracle = NetworkHistogram.oracle_estimator(
    A, oracle_latents, NetworkHistogram.UnitIntervalConvertor(n_bins));

starting_labels = shuffle(oracle_latents);

max_iter = 1_000_000
stalled_iters = 10_000

res_new = NetworkHistogram.nethist_continuous(
    A, k,
    starting_labels;
    bins = n_bins
);

ENV["JULIA_CONDAPKG_VERBOSITY"] = "-1" # hide conda messages #hide
using PythonCall

θ_oracle = Graphons._extract_param.(res_oracle.model.θ);
θ_hat = Graphons._extract_param.(res_new.model.θ);
perm, plan = NetworkHistogram.get_perm_alignment(θ_oracle, θ_hat);

let
    fig = Mke.Figure(size = (600, 400))
    ax = Mke.Axis(fig[1, 1], title = "OT plan heatmap",
        xlabel = "Fitted groups", ylabel = "Oracle groups")
    hm = Mke.heatmap!(ax, plan, colormap = :binary)
    Mke.Colorbar(fig[1, 2], hm)
    fig
    Mke.display(fig) #src
end

fitted_labels = map(x -> perm[x], res_new.labels);
res_ot_aligned = NetworkHistogram.oracle_estimator(
    A, fitted_labels, NetworkHistogram.UnitIntervalConvertor(n_bins),
    name = "aligned with OT perm");
# NetworkHistogram.align_res_true_latents!(res_new, res_oracle.labels);

xs = range(0, 1; length = 20)

function viz_one_group!(axis, g1, g2, A, ξs, res_oracle, res_new, xs; n_viz = 20, p = p)
    nodes_1 = findall(res_oracle.labels .== g1)
    nodes_2 = findall(res_oracle.labels .== g2)
    edge_values = [A[x, y] for y in nodes_2 for x in nodes_1]
    Mke.vlines!(axis, edge_values, ymax = 0.025, color = :lightgray)
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

fig = Mke.Figure(size = (1000, 1000))
for g in 1:k
    for g2 in 1:g
        ax = Mke.Axis(fig[g, g2])#, title = "Group $g vs Group $g2", xlabel = "Edge Value",ylabel = "Density")
        Mke.hidedecorations!(ax)
        viz_one_group!(ax, g, g2, A, ξs, res_oracle,
            res_ot_aligned, xs, p = p, n_viz = 5)
    end
end
fig
Mke.display(fig) #src

##

clustering_res = kmeans(A, k)

res_kmeans = NetworkHistogram.oracle_estimator(
    A, assignments(clustering_res), NetworkHistogram.UnitIntervalConvertor(n_bins);
    type_suff_stats = Val(:categorical),
    name = "k-means");

NetworkHistogram.align_res_true_latents!(res_kmeans, res_oracle.labels);

fig = Mke.Figure(size = (1000, 1000))
for g in 1:k
    for g2 in 1:g
        ax = Mke.Axis(fig[g, g2], title = "Group $g vs Group $g2", xlabel = "Edge Value",
            ylabel = "Density")
        viz_one_group!(
            ax, g, g2, A, ξs, res_oracle, res_kmeans, xs, p = p, n_viz = 5)
    end
end

fig
display(fig) #src
