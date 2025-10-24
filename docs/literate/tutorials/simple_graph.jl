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
n = 2000
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

# First, we need to represent our graph in a format that the package understands.
# We can use an `EdgeList` to store the edges of the graph.
edge_list = EdgeList(A);

# We also need to define the model for the edges. Since our graph is unweighted,
# we can use a `Bernoulli` distribution. The `Dist` wrapper is used to
# handle aggregation of distributions.
import NetworkHistogram: Dist, Assignment, nethist
dist = NetworkHistogram.Bernoulli(0.5) # The initial probability doesn't matter much.

# We start with a random initial assignment of nodes to `k=5` groups.
k = 10
oracle_labels = ordered_start_labels(n, k);

initial_assignment = shuffle(oracle_labels);

# Now, we create an `Assignment` object, which holds all the information
# about the model and the current state of the node groupings.
oracle_estimator = Assignment(oracle_labels, edge_list, Dist(dist));
sbm_oracle = NetworkHistogram.to_block_model(oracle_estimator);
Mke.heatmap(sbm_oracle, colormap = :binary, colorrange = (0, 1))

println("Log-likelihood of oracle estimator: ", loglikelihood(oracle_estimator))
# `NetworkHistogram.jl` provides optimization algorithms to improve the initial assignment.
# Let's use the `nethist` function with `GreedyParams`, which iteratively moves nodes between
# groups to maximize the log-likelihood.

# params_opti = NetworkHistogram.GreedyParams(
#     100_000, NetworkHistogram.RandomNodeSwap(), NetworkHistogram.Strict(),
#     NetworkHistogram.PreviousBestValue(2_000), false);

# a = nethist(A, dist, initial_assignment, params_opti, false);

res = NetworkHistogram.nethist_binary_edges(A,
    initial_assignment, GreedyParams(
        1_000_000,
        RandomGroupSwap(),
        Strict(),
        PreviousBestValue(1_000, Inf, :min),
        true
    ));

a = Assignment(res.node_labels, edge_list, Dist(dist));
println("Log-likelihood after optimization: ", loglikelihood(a))

# The `Assignment` object `a` now contains the optimized node groupings and
# the fitted network histogram parameters.

# We can visualize the fitted histogram.
heatmap_params(a, ordering = false, colorrange = (0, 1))

# We can convert it to a block model for easier interpretation.

# res = NethistResult(a);

let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["True Graphon W(u,v)", "Oracle Estimator", "Fitted Network Histogram"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], w, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], sbm_oracle,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], res.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

# the block labels found by the optimization are not necessarily aligned with the true latent positions, hence the need to align them for better visualization.

NetworkHistogram.align_res_true_latents!(res, a, oracle_estimator.node_labels);

# and display the true function, the oracle estimator, and the fitted model
let
    fig = Mke.Figure(size = (1220, 400))
    titles = ["True Graphon W(u,v)", "Oracle Estimator", "Fitted Network Histogram"]
    axes = [Mke.Axis(fig[1, i], aspect = Mke.DataAspect(), title = titles[i]) for i in 1:3]
    Mke.heatmap!(axes[1], w, colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[2], sbm_oracle,
        colormap = :binary, colorrange = (0, 1))
    Mke.heatmap!(axes[3], res.model, colormap = :binary, colorrange = (0, 1))
    Mke.Colorbar(fig[1, 4], colormap = :binary,
        limits = (0, 1), label = "Edge Probability", width = 20)
    fig
end

# We can even fit a Stochastic Shape Model quite easily from the fitted SBM.

using Clustering

ξ = NetworkHistogram.node_labels_to_latents(res.node_labels, res.model);
shape_range = 1:(k * (k + 1) ÷ 2 - 1)
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res.model, A, ξ, shape_range)

using Kneedle
kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing)
#  Let's extract the optimal number of shapes using the Kneedle algorithm:

k_knee = knees(kr)[1]
ssm_knee = SSM(res.model, k_knee)

Mke.heatmap(ssm_estimated, colormap = :binary, colorrange = (0, 1))
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
