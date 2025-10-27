import NetworkHistogram: SuffStats, add_sample, remove_sample, make_k_block, score
using StaticArrays
using Accessors

struct MyCustomSuffStats{M, T} <: SuffStats
    h::SVector{M, T}
end

function MyCustomSuffStats(num_categories::Int)
    h = SVector{num_categories, Int}(zeros(Int, num_categories))
    return MyCustomSuffStats{num_categories, Int}(h)
end

@inline function add_sample(ss::MyCustomSuffStats, sample::Int)
    ss = @set ss.h[sample] += 1
    return ss
end

@inline function remove_sample(ss::MyCustomSuffStats, sample::Int)
    ss = @set ss.h[sample] -= 1
    return ss
end

function make_k_block(k, ::Val{:custom}; num_categories, kwargs...)
    k_block = SymArray{MyCustomSuffStats{num_categories, Int}}(undef, k, k)
    fill!(k_block, MyCustomSuffStats(num_categories))
    return k_block
end

@inline function score(ss::MyCustomSuffStats; kwargs...)
    n = sum(ss.h)
    return n - sum(abs2, ss.h) / max(n, 1)
end

##

using Distributions
using NetworkHistogram
using Random
using StatsBase

function W_multiplex(x, y)
    ps = zeros(4)
    ps[2] = sqrt(abs(x - y)) / 2           # layer 1 only
    ps[3] = abs(sin(2π * x) * sin(2π * y)) / 2  # layer 2 only
    ps[4] = min(x, y) / 2                   # both layers
    ps[1] = 1 - sum(ps[2:4])                # no edge
    return DiscreteNonParametric(0:3, SVector{4}(ps))
end

m = 4
graphon = DecoratedGraphon(W_multiplex)

n = 2000
true_latents = range(0, 1; length = n)
A = sample_graph(graphon, true_latents);

k = 20
oracle_labels = ordered_start_labels(n, k);
initial_labels = shuffle(oracle_labels);

max_iter = 1_000_000
stalled_iters = 5_000

data = A .+ 1;  # shift to 1,2,3,4 for categorical
es_new = NetworkHistogram.GreedySuffStats(data, initial_labels, num_categories = m,
    type_suff_stats = :custom,
    max_iter = max_iter,
    swap_rule = NetworkHistogram.RandomGroupSwap(),
    stop_rule = NetworkHistogram.PreviousBestValue(stalled_iters, Inf, :min),
    progress = true,
    dist = Categorical(m)
);
node_labels_es_new = NetworkHistogram.estimate!(
    es_new, data, initial_labels; dist = Categorical(m),
    iter_progress = 10_000)

function params(ss::Union{NetworkHistogram.CategoricalSuffStats, MyCustomSuffStats})
    ss.h ./ sum(ss.h)
end
parameters = Matrix{SVector{m, Float64}}(undef, k, k)
@inbounds for j in 1:k, i in 1:k
    parameters[i, j] = SVector{m, Float64}(
        params(es_new.block_ss[i, j])...)
end
model_es_new = NetworkHistogram.DecoratedSBM(
    DiscreteNonParametric.(Ref(0:(m - 1)), parameters), counts(node_labels_es_new) ./
                                                        length(node_labels_es_new));

res_new = NetworkHistogram.NethistResult(node_labels_es_new, model_es_new);
NetworkHistogram.align_res_true_latents!(res_new, oracle_labels);
