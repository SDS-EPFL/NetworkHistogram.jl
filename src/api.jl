function nethist_categorical(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        params::GreedyParams = GreedyParams())
    convertor = CategoricalConvertor(A)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    _nethist(
        A, labels_start,
        convertor,
        Val(:categorical),
        params,
        num_categories = num_bins(convertor)
    )
end

function nethist_continuous(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        num_bins_::Int = 10,
        params::GreedyParams = GreedyParams())
    convertor = UnitIntervalConvertor(num_bins_)
    @info "Using $(num_bins(convertor)) discrete categories for edge values"
    _nethist(
        A, labels_start,
        convertor,
        Val(:categorical),
        params,
        num_categories = num_bins(convertor)
    )
end

function nethist_binary(
        A, k,
        labels_start = ordered_start_labels(size(A, 1), k);
        params::GreedyParams = GreedyParams())
    _nethist(
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
    es = GreedySuffStats(
        data, labels_start,
        type_suff_stats = type_suff_stats,
        max_iter = params.max_iter,
        swap_rule = params.node_swap_rule,
        stop_rule = params.stop_rule,
        progress = params.display_progress;
        kwargs...
    )
    node_labels, parameters = estimate!(
        es, data, labels_start; iter_progress = params.progress_freq)

    return convert_to_result(node_labels, convertor, parameters)
end

function oracle_estimator(
        data, oracle_labels, convertor; type_suff_stats = Val(:categorical))
    k = length(unique(oracle_labels))
    # allocate sufficient statistics blocks
    block_ss = make_k_block(k, type_suff_stats; data = data)
    block_ss_swap = make_k_block(k, type_suff_stats; data = data)
    es_dummy = GreedySuffStats(block_ss, block_ss_swap, RandomGroupSwap(),
        PreviousBestValue(1_000, Inf, :min), 1)
    init!(es_dummy, data, oracle_labels)
    @info "Oracle estimator loss: $(loss(es_dummy, norm = get_num_obs(data)))"
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
