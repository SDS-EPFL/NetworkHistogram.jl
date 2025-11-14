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
    ps[4] = min(x, y) / 2                   # both layers
    ps[1] = 1 - sum(ps[2:4])                # no edge
    return DiscreteNonParametric(0:3, SVector{4}(ps))
end

function W3(x, y)
    ps = zeros(4)
    ps[1] = 3 * x * y
    ps[2] = 3 * sin(2 * π * x) * sin(2 * π * y)
    ps[3] = exp(-3 * (x - 0.5)^2 - 3 * (y - 0.5)^2)
    ps[4] = 2 - 3 * (x + y)
    e_ps = exp.(ps)
    return DiscreteNonParametric(0:3, SVector{4}(e_ps ./ sum(e_ps)))
end

graphon = DecoratedGraphon(W3)

let
    fig = Mke.Figure(size = (4 * h, h))
    for m in 1:4
        ax = Mke.Axis(fig[1, m], aspect = Mke.DataAspect())
        Mke.heatmap!(ax, graphon, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
    display(fig) #src
end

n = 1000
true_latents = range(0, 1; length = n)
A = sample_graph(graphon, true_latents);

k = 14
oracle_labels = ordered_start_labels(n, k);
initial_labels = shuffle(oracle_labels);

oracle_res = NetworkHistogram.oracle_estimator(
    A, oracle_labels, NetworkHistogram.CategoricalConvertor(A));

res = NetworkHistogram.nethist_categorical(A, k, initial_labels)

# Visualize the fitted models for different numbers of groups after aligning with true latents

NetworkHistogram.align_res_true_latents!(res, oracle_res.labels);
let
    fig = Mke.Figure(size = (4 * h, h))
    for m in 1:4
        ax = Mke.Axis(fig[1, m], aspect = Mke.DataAspect())
        Mke.heatmap!(ax, res.model, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
    display(fig) #src
end

# We can also align the fitted model to the true one using optimal transport. We need to load the `PythonCall.jl`
# package for that, as we will use the `POT` Python library.

ENV["JULIA_CONDAPKG_VERBOSITY"] = "-1" # hide conda messages #hide
using PythonCall
θ_oracle = probs.(oracle_res.model.θ);
θ_hat = probs.(res.model.θ);

perm = NetworkHistogram.get_perm_alignment(θ_hat, θ_oracle);

θ_hat_aligned = θ_hat[perm, perm];
estimator_aligned = DecoratedSBM(DiscreteNonParametric.(Ref(0:3), θ_hat_aligned),
    res.model.size[perm]);

let
    fig = Mke.Figure(size = (2 * h, h))
    for m in 1:4
        ax = Mke.Axis(
            fig[1, m], aspect = Mke.DataAspect(), ylabel = m == 1 ? "Estimated" : "")
        Mke.heatmap!(ax, estimator_aligned, k = m, colormap = :binary, colorrange = (0, 1))
        ax2 = Mke.Axis(
            fig[2, m], aspect = Mke.DataAspect(), ylabel = m == 1 ? "Oracle" : "")
        Mke.heatmap!(
            ax2, oracle_res.model, k = m, colormap = :binary, colorrange = (0, 1))
    end
    fig
    display(fig) #src
end

# The fitted network histogram can be further processed to obtain a smoother estimate of the underlying graphon.

using Clustering
shape_range = 1:30
ssm_estimated, criterion_values = Graphons.estimate_ssm(
    res.model, A, true_latents, shape_range);

using Kneedle
kr = kneedle(shape_range, criterion_values, "convex_dec", 1, scan_type = :smoothing);
#  Let's extract the optimal number of shapes using the Kneedle algorithm:

k_knee = knees(kr)[1]
ssm = SSM(res.model, k_knee)

models_to_plot = [graphon, res.model, ssm_estimated, ssm]
model_names = ["True graphon", "Block model",
    "SSM argmin k=$(length(ssm_estimated.θ))", "SSM knee k=$k_knee"]

let
    fig = Mke.Figure(size = (4 * h, length(models_to_plot) * h))
    for (i, model) in enumerate(models_to_plot)
        for m in 1:4
            ax = Mke.Axis(
                fig[i, m], aspect = Mke.DataAspect(), ylabel = m == 1 ? model_names[i] : "")
            Mke.hidedecorations!(ax, label = false)
            Mke.heatmap!(ax, model, k = m, colormap = :lipari, colorrange = (0, 1))
        end
    end
    Mke.Colorbar(fig[2:3, end + 1], colormap = :lipari,
        limits = (0, 1), width = 0.05 * h)
    fig
    display(fig) #src
end
