#=
#  A Simple Graphon Tutorial with NetworkHistogram.jl
=#

# This tutorial introduces the concept of a graphon, demonstrates how to sample a graph from one, and then shows how to estimate the graphon from the sampled graph using the Network Histogram method provided by `NetworkHistogram.jl`.

# ## What is a Graphon?

# A graphon (or graph function) is a symmetric, measurable function $$W: [0, 1]^2 \to [0, 1]$$.

# It serves as a generative model for random graphs. Think of it as a continuous and more general version of a stochastic block model.

# In simple terms, each node `i` in a graph is assigned a latent (unobserved) position $u_i \in [0, 1]$. The probability of an edge existing between two nodes `i` and `j` is then given by the graphon function evaluated at their latent positions:

# Let's define a simple graphon in Julia. For this example, we'll use a step-function-like graphon that resembles a stochastic block model.

import CairoMakie as Mke
using LinearAlgebra
using Random
import StatsBase: inverse_rle
using Statistics
using NetworkHistogram
using Distributions

h = 300; # hide
Random.seed!(1234);

# Define a simple step-function graphon
w = SimpleContinuousGraphon((x, y) -> x * y)

# We can visualize this graphon as a heatmap.
let
    fig = Mke.Figure(size = (h + 20, h))
    ax = Mke.Axis(fig[1, 1], title = "True Graphon W(u,v)")
    hm = Mke.heatmap!(ax, w, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 2], hm)
    fig
end

#md
# ## Sampling a Graph from a Graphon

# To generate a random graph from a graphon, we follow these steps:
# 1.  **Assign latent positions:** For a graph with `n` nodes, we sample `n` independent and identically distributed random variables $u_1, u_2, \dots, u_n$ from a Uniform(0, 1) distribution. These are the latent positions of our nodes.
# 2.  **Generate edges:** For each pair of nodes `(i, j)` with `i < j`, we generate a random number from a Bernoulli distribution with probability $W(u_i, u_j)$. This determines whether an edge exists between them. The resulting adjacency matrix `A` will be symmetric.
# Let's sample a graph with 2000 nodes from our graphon `W`.
n = 3000
u_true = rand(n);  # Latent positions
A = sample_graph(w, u_true);

# We can visualize the adjacency matrix of the sampled graph.
# To make the block structure visible, we sort the nodes by their latent positions.
perm = sortperm(u_true)
A = A[perm, perm]
let
    fig = Mke.Figure(size = (h, h))
    ax = Mke.Axis(
        fig[1, 1], title = "Sampled Adjacency Matrix (Sorted)", aspect = Mke.DataAspect())
    Mke.heatmap!(ax, A, colormap = :binary)
    fig
end

#md
# ## The Network Histogram Method

# The Network Histogram method is a non-parametric approach to estimate a graphon from a single observed network. The core idea is to approximate the (unknown) graphon `W` with a piecewise constant function.

# This is achieved by:
# 1.  **Partitioning the nodes:** The nodes of the graph are partitioned into `k` groups.
# 2.  **Estimating block probabilities:** The probability of an edge between any two groups is estimated by the density of edges between them.
# 3.  **Constructing the histogram:** These estimated probabilities form a `k x k` matrix, which is a step-function approximation of the true graphon.

# The main challenge is to find the optimal partition of nodes. `NetworkHistogram.jl` provides tools to find a good partition by optimizing an objective function, such as the log-likelihood of the observed graph under the model.

# ## Fitting a Network Histogram with NetworkHistogram.jl

# Now, let's use `NetworkHistogram.jl` to fit a network histogram to the graph `A` we sampled earlier. We will try to recover the underlying 2-block structure.

# We start with a random initial assignment of nodes to `k=5` groups.
k = 10
oracle_labels = ordered_start_labels(n, k);

initial_assignment = shuffle(oracle_labels);

## We can compute the "oracle" estimator, which uses the true latent positions to assign nodes to groups. This serves as a benchmark for our estimation.
oracle_res = NetworkHistogram.oracle_estimator(
    A, oracle_labels, NetworkHistogram.BinaryConvertor(); type_suff_stats = Val(:binary));

let
    fig = Mke.Figure(size = (400, 300))
    ax = Mke.Axis(fig[1, 1], aspect = Mke.DataAspect())
    Mke.heatmap!(ax, oracle_res.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 2], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end
##
# `NetworkHistogram.jl` provides optimization algorithms to improve the initial assignment.
# Let's use the `nethist` function with `GreedyParams`, which iteratively moves nodes between
# groups to maximize the log-likelihood.

# params_opti = NetworkHistogram.GreedyParams(
#     100_000, NetworkHistogram.RandomNodeSwap(), NetworkHistogram.Strict(),
#     NetworkHistogram.PreviousBestValue(2_000), false);

# a = nethist(A, dist, initial_assignment, params_opti, false);

res = NetworkHistogram.nethist_binary(A, k, initial_assignment);

let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["True Graphon W(u,v)", "Oracle Estimator", "Fitted Network Histogram"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], w, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], oracle_res.model,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], res.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

# the block labels found by the optimization are not necessarily aligned with the true latent positions, hence the need to align them for better visualization.

NetworkHistogram.align_res_true_latents!(res, oracle_res.labels);

# and display the true function, the oracle estimator, and the fitted model
let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["True Graphon W(u,v)", "Oracle Estimator", "Fitted Network Histogram"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], w, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], oracle_res.model,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], res.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

# We can even fit a Stochastic Shape Model quite easily from the fitted SBM.

using Clustering

# ξ = NetworkHistogram.node_labels_to_latents(res.labels, res.model);
shape_range = 1:(k * (k + 1) ÷ 2 - 1)
ssm_estimated,
criterion_values = Graphons.estimate_ssm(
    res.model, A, res.labels, shape_range)

using Kneedle
kr = kneedle(shape_range, criterion_values, "convex_dec", 1,
    kneedle_scan_algorithm = ScanSmoothing(; S = 1.0))
#  Let's extract the optimal number of shapes using the Kneedle algorithm:

k_knee = knees(kr)[1]
ssm_knee = SSM(res.model, k_knee)

println("Number of shapes in SSM argmin: ", length(ssm_estimated.θ))
println("Number of shapes in SSM knee: ", length(ssm_knee.θ))
println("Number of shapes in SBM: ", length(res.model.θ))

# We greatly reduced the number of parameters from the original SBM estimate while preserving much of the structure of the estimated graphon as seen below:

let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["SBM", "SSM argmin", "SSM knee"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], res.model, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], ssm_estimated,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], ssm_knee, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

##

k_kmeans = 10;
clustering_res = kmeans(A, k_kmeans);

res_kmeans = NetworkHistogram.oracle_estimator(
    A, assignments(clustering_res), NetworkHistogram.BinaryConvertor();
    type_suff_stats = Val(:binary),
    name = "k-means");

NetworkHistogram.align_res_true_latents!(res_kmeans, oracle_res.labels);

# and display the true function, the oracle estimator, and the fitted model
let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["True Graphon W(u,v)", "Oracle Estimator", "Fitted Network Histogram"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], w, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], oracle_res.model,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], res_kmeans.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

# ξ = NetworkHistogram.node_labels_to_latents(res.labels, res.model);
shape_range = 1:(k_kmeans * (k_kmeans + 1) ÷ 2 - 1)
ssm_estimated,
criterion_values = Graphons.estimate_ssm(
    res_kmeans.model, A, res_kmeans.labels, shape_range)

using Kneedle
kr = kneedle(shape_range, criterion_values, "convex_dec", 1,
    kneedle_scan_algorithm = ScanSmoothing(; S = 1.0))
#  Let's extract the optimal number of shapes using the Kneedle algorithm:

k_knee = knees(kr)[1]
ssm_knee = SSM(res_kmeans.model, k_knee)

println("Number of shapes in SSM argmin: ", length(ssm_estimated.θ))
println("Number of shapes in SSM knee: ", length(ssm_knee.θ))
println("Number of shapes in SBM: ", length(res_kmeans.model.θ))

# We greatly reduced the number of parameters from the original SBM estimate while preserving much of the structure of the estimated graphon as seen below:

let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["SBM", "SSM argmin", "SSM knee"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], res_kmeans.model, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], ssm_estimated,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], ssm_knee, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end
