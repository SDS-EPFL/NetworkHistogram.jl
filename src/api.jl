
function nethist_categorical(
        A, k,
        labels_start = shuffle(ordered_start_labels(size(A, 1), k));
        max_iter = 1_000_000,
        stalled_iters = 5_000,
        progress_bar = true)
    _nethist(
        A, labels_start,
        CategoricalConvertor(A),
        :categorical;
        max_iter = max_iter,
        stalled_iters = stalled_iters,
        progress_bar = progress_bar
    )
end

function nethist_continuous(
        A, k,
        labels_start = shuffle(ordered_start_labels(size(A, 1), k));
        num_bins_::Int = 10,
        max_iter = 1_000_000,
        stalled_iters = 5_000,
        progress_bar = true)
    _nethist(
        A, labels_start,
        UnitIntervalConvertor(num_bins_),
        :categorical;
        max_iter = max_iter,
        stalled_iters = stalled_iters,
        progress_bar = progress_bar
    )
end

function _nethist(A, labels_start, convertor, type_suff_stats; max_iter = max_iter,
        stalled_iters = 10_000,
        progress_bar = true)
    data = convertor.(A)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    es = NetworkHistogram.GreedySuffStats(
        data, labels_start, num_categories = num_bins(convertor),
        type_suff_stats = type_suff_stats,
        max_iter = max_iter,
        swap_rule = NetworkHistogram.RandomGroupSwap(),
        stop_rule = NetworkHistogram.PreviousBestValue(stalled_iters, Inf, :min),
        progress = progress_bar
    )
    node_labels, parameters = NetworkHistogram.estimate!(
        es, data, labels_start; iter_progress = 10_000)
    model = NetworkHistogram.DecoratedSBM(to_distribution.(convertor, parameters),
        counts(node_labels) ./ length(node_labels))
    return NetworkHistogram.NethistResult(node_labels, model)
end

# function nethist_binary_edges(A, initial_node_labels, params::GreedyParams)
#     k = length(unique(initial_node_labels))
#     data, counts_main, counts_swap, realized, realized_swap = prepare_data_cat(A, k)
#     es = GreedyAverage(
#         counts_main, counts_swap, realized, realized_swap,
#         params.max_iter, params.swap_rule, params.stop_rule)
#     node_labels = estimate(es, data, initial_node_labels, progress = params.progress_bar)
#     sizes = counts(node_labels) ./ length(node_labels)

#     θ = Matrix{Float64}(undef, k, k)
#     @inbounds for j in 1:k, i in 1:k
#         θ[i, j] = es.realized[i, j][2] / max(1, es.counts[i, j])
#     end
#     model = SBM(θ, sizes)

#     return NethistResult(node_labels, model)
# end

# function nethist_discrete_edges(
#         A, initial_node_labels, params::GreedyParams, m = length(unique(A)))
#     k = length(unique(initial_node_labels))
#     data, counts_main, counts_swap, realized, realized_swap = prepare_data_cat(A, k, m = m)
#     es = GreedyAverage(
#         counts_main, counts_swap, realized, realized_swap,
#         params.max_iter, params.swap_rule, params.stop_rule)
#     node_labels = estimate(es, data, initial_node_labels, progress = params.progress_bar)
#     sizes = counts(node_labels) ./ length(node_labels)
#     parameters = Matrix{SVector{m, Float64}}(undef, k, k)
#     @inbounds for j in 1:k, i in 1:k
#         parameters[i, j] = [es.realized[i, j][c] / max(es.counts[i, j], 1) for c in 1:m]
#     end
#     s = zero(eltype(A)) in A ? collect(0:(m - 1)) : 1:m
#     model = DecoratedSBM(DiscreteNonParametric.(Ref(s), parameters), sizes)
#     return NethistResult(node_labels, model), es
# end

# function nethist_continuous_edges(A_cont, initial_node_labels, params::GreedyParams;
#         num_bins_::Int = 10, lower_bound = quantile(A_cont[:], 0.01), upper_bound = quantile(
#             A_cont[:], 0.99))
#     convertor = ContinuousConvertor(lower_bound, upper_bound, num_bins_)
#     A = convertor.(A_cont)
#     @info "Discretized continuous edge values into $(num_bins(convertor)) bins"
#     res_cat = nethist_discrete_edges(
#         A, initial_node_labels, params, num_bins(convertor))
#     parameters = NetworkHistogram.HistDistribution.(
#         Graphons._extract_param.(res_cat.model.θ), convertor)
#     model = DecoratedSBM(parameters, res_cat.model.size)
#     return NethistResult(res_cat.node_labels, model), res_cat, A
# end
