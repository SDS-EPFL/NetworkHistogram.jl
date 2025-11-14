function nethist_categorical(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        params::GreedyParams = GreedyParams(stalled_iters = 50_000, max_iter = 2_000_000))
    convertor = CategoricalConvertor(A)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    return _nethist(
        A, labels_start,
        convertor,
        Val(:categorical),
        params;
        num_categories = num_bins(convertor)
    )
end

function nethist_continuous(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        bins::Int = 10,
        params::GreedyParams = GreedyParams())
    convertor = UnitIntervalConvertor(bins)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    return _nethist(
        A, labels_start,
        convertor,
        Val(:categorical),
        params;
        num_categories = num_bins(convertor)
    )
end

function nethist_binary(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        params::GreedyParams = GreedyParams())
    return _nethist(
        A, labels_start,
        BinaryConvertor(),
        Val(:binary),
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
    es = make_greedy_suffstats_estimator(
        data,
        labels_start;
        type_suff_stats = type_suff_stats,
        max_iter = params.max_iter,
        node_swap_rule = params.node_swap_rule,
        stop_rule = params.stop_rule,
        kwargs...
    )
    node_labels, parameters = estimate!(
        es, data, labels_start;
        progress = params.display_progress,
        iter_progress = params.progress_freq)

    return convert_to_result(node_labels, convertor, parameters)
end

function oracle_estimator(
        A, oracle_labels, convertor; type_suff_stats = Val(:categorical),
        name = "oracle", kwargs...)

    # prepare data
    k = length(unique(oracle_labels))
    data = convertor.(A)
    # prepare suff stats
    block_ss = make_k_block(
        k, type_suff_stats; data = data, num_categories = num_bins(convertor), kwargs...)

    # compute oracle suff stats
    es_dummy = GreedySuffStats(block_ss, copy(block_ss), RandomGroupSwap(),
        PreviousBestValue(1_000, Inf, :min), 1)
    init!(es_dummy, data, oracle_labels)

    # retrieve parameters
    @info "$name estimator loss: $(loss(es_dummy, norm = get_num_obs(data)))"
    parameters = to_params.(es_dummy.block_ss)
    return convert_to_result(oracle_labels, convertor, parameters)
end

function convert_to_result(node_labels, convertor, parameters)
    model = DecoratedSBM(to_distribution.(convertor, parameters),
        counts(node_labels) ./ length(node_labels))
    return NethistResult(node_labels, model)
end

function convert_to_result(node_labels, convertor::BinaryConvertor, parameters)
    model = SBM(to_distribution.(convertor, parameters),
        counts(node_labels) ./ length(node_labels))
    return NethistResult(node_labels, model)
end
