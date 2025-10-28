function nethist_categorical(
        A, k,
        labels_start = shuffle(ordered_start_labels(size(A, 1), k));
        params::GreedyParams = GreedyParams())
    _nethist(
        A, labels_start,
        CategoricalConvertor(A),
        Val(:categorical),
        params
    )
end

function nethist_continuous(
        A, k,
        labels_start = shuffle(ordered_start_labels(size(A, 1), k));
        num_bins_::Int = 10,
        params::GreedyParams = GreedyParams())
    _nethist(
        A, labels_start,
        UnitIntervalConvertor(num_bins_),
        Val(:categorical),
        params
    )
end

function _nethist(
        A, labels_start, convertor, type_suff_stats,
        params::GreedyParams = GreedyParams(); kwargs...)
    if !params.warm_start
        reset!(params)
    end
    data = convertor.(A)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    es = NetworkHistogram.GreedySuffStats(
        data, labels_start, num_categories = num_bins(convertor),
        type_suff_stats = type_suff_stats,
        max_iter = params.max_iter,
        swap_rule = params.node_swap_rule,
        stop_rule = params.stop_rule,
        progress = params.display_progress;
        kwargs...
    )
    node_labels, parameters = NetworkHistogram.estimate!(
        es, data, labels_start; iter_progress = params.progress_freq)
    model = NetworkHistogram.DecoratedSBM(to_distribution.(convertor, parameters),
        counts(node_labels) ./ length(node_labels))
    return NetworkHistogram.NethistResult(node_labels, model)
end
