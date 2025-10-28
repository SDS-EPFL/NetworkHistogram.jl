abstract type SBMEstimator end

abstract type Result end

struct NethistResult{L, M} <: Result
    labels::L
    model::M
end

struct GreedySuffStats{M, NodeR <: NodeSwapRule, StopR <: StopRule} <: SBMEstimator
    block_ss::M
    block_ss_swap::M
    node_swap_rule::NodeR
    stop_rule::StopR
    max_iter::Int
end

function init!(es::GreedySuffStats, data, node_labels)
    # Initialize the sufficient statistics for each block
    @inbounds for j in axes(data, 2)
        gj = node_labels[j]
        for i in 1:(j - 1)  # More efficient than i < j check inside loop
            edge_value = data[i, j]
            gi = node_labels[i]
            es.block_ss[gi, gj] = add_sample(es.block_ss[gi, gj], edge_value, i, j)
            es.block_ss_swap[gi, gj] = add_sample(
                es.block_ss_swap[gi, gj], edge_value, i, j)
        end
    end
end

# TODO: allow for non-symmetric data
@inline function loss(matrix_ss::SymArray{<:SuffStats}; norm = 1.0)
    total_loss = 0.0
    for m in matrix_ss.uppertrian.nzval
        total_loss += loss(m)
    end
    return total_loss / norm
end

@inline function loss(matrix_ss::AbstractMatrix{<:SuffStats}; norm = 1.0)
    total_loss = 0.0
    @inbounds for j in axes(matrix_ss, 2)
        for i in 1:j
            inter = loss(matrix_ss[i, j])
            total_loss += inter
        end
    end
    return total_loss / norm
end

function GreedySuffStats(
        data, node_labels; type_suff_stats = Val(:categorical), max_iter = 10000,
        node_swap_rule = RandomGroupSwap(), stop_rule = PreviousBestValue(5_000, Inf, :min),
        dist = nothing,
        kwargs...)
    # derive user input
    k = length(unique(node_labels))

    # allocate sufficient statistics blocks
    block_ss = make_k_block(k, type_suff_stats; data = data, dist = dist, kwargs...)
    block_ss_swap = make_k_block(k, type_suff_stats; data = data, dist = dist, kwargs...)

    # create estimator
    return GreedySuffStats{typeof(block_ss), typeof(node_swap_rule), typeof(stop_rule)}(
        block_ss, block_ss_swap, node_swap_rule, stop_rule, max_iter)
    return es
end

function estimate!(
        es::GreedySuffStats,
        data,
        node_labels_init;
        progress = true,
        iter_progress = 5000
)
    # Initialize node labels
    node_labels = copy(node_labels_init)
    n = length(node_labels)
    k = length(unique(node_labels))
    n_edges = n * (n - 1) / 2
    init!(es, data, node_labels)

    # Progress tracking
    pbar = ProgressUnknown(
        enabled = progress,
        showspeed = true,
        desc = "Greedy search: "
    )

    # Update progress bar only every N iterations to reduce overhead
    progress_update_interval = max(1, es.max_iter ÷ iter_progress)
    # Initial log-likelihood

    current_loss = loss(es.block_ss, norm = n_edges)
    reset!(es.stop_rule, current_loss)
    # Main optimization loop
    for iter in 1:(es.max_iter)
        # Select two nodes to potentially swap
        index1, index2 = select_indices_swap(node_labels, es.node_swap_rule, k)

        group1 = node_labels[index1]
        group2 = node_labels[index2]

        @inbounds for j in axes(data, 2)
            if j != index1 && j != index2
                # extract data
                groupj = node_labels[j]
                edge_value_1 = data[j, index1]
                edge_value_2 = data[j, index2]

                es.block_ss_swap[group1, groupj] = remove_sample(
                    es.block_ss_swap[group1, groupj], edge_value_1, j, index1)
                es.block_ss_swap[group2, groupj] = add_sample(
                    es.block_ss_swap[group2, groupj], edge_value_1, j, index1)

                es.block_ss_swap[group2, groupj] = remove_sample(
                    es.block_ss_swap[group2, groupj], edge_value_2, j, index2)
                es.block_ss_swap[group1, groupj] = add_sample(
                    es.block_ss_swap[group1, groupj], edge_value_2, j, index2)
            end
        end

        # tentative swap
        @inbounds node_labels[index1], node_labels[index2] = group2, group1
        new_loss = loss(es.block_ss_swap, norm = n_edges)

        if compare_to_best(new_loss, current_loss, es.stop_rule)
            # apply swap
            copy!(es.block_ss, es.block_ss_swap)
            current_loss = new_loss
        else
            # revert swap
            node_labels[index1], node_labels[index2] = group1, group2
            # revert sufficient statistics
            copy!(es.block_ss_swap, es.block_ss)
        end

        if progress && (iter % progress_update_interval == 0 || iter == es.max_iter)
            update!(
                pbar, iter;
                showvalues = [
                    ("loss", current_loss),
                    info_to_print(es.stop_rule)
                ])
        end

        # Check stopping criterion
        if stopping_rule(current_loss, es.stop_rule)
            break
        end
    end
    finish!(pbar)
    @info "Optimization finished. Final loss: $current_loss"

    return node_labels, to_params.(es.block_ss)
end
