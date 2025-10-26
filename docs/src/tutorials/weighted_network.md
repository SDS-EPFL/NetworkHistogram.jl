```@meta
EditURL = "../../literate/tutorials/weighted_network.jl"
```

# Decorated Graphon Tutorial for Weighted Networks

````@example weighted_network
using NetworkHistogram
using Distributions
import CairoMakie as Mke
using LinearAlgebra

graphon = DecoratedGraphon((x, y) -> Exponential(3 * x * y + 1))

let
    fig = Mke.Figure(size = (600, 300))
    ax = Mke.Axis(fig[1, 1], aspect = Mke.DataAspect())
    Mke.heatmap!(ax, graphon, colormap = :viridis)
    fig
end

n = 500
k = 10
A = sample_graph(graphon, n) .* Symmetric(rand(Bernoulli(0.7), n, n));
oracle_latents = ordered_start_labels(n, k);
starting_labels = shuffle(oracle_latents);

res = NetworkHistogram.nethist_continuous_edges(A,
    starting_labels, GreedyParams(
        1_000_000,
        RandomGroupSwap(),
        Strict(),
        PreviousBestValue(5_000, Inf, :min),
        true # progress bar
    ));
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

