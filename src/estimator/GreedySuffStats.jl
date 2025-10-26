abstract type SuffStats end

function add_sample(ss::SuffStats, sample)
    @error("add_sample not implemented for $(typeof(ss)) and sample $(typeof(sample)) \n
        you may need to implement a custom sufficient statistics type")
end
function remove_sample(ss::SuffStats, sample)
    @error("remove_sample not implemented for $(typeof(ss)) and sample $(typeof(sample)) \n
        you may need to implement a custom sufficient statistics type")
end

function make_k_block(k, suff_stats_type; kwargs...)
    @error("make_k_block not implemented for sufficient statistics type $(suff_stats_type) \n
        you may need to implement a custom sufficient statistics type")
end

function score(ss::SuffStats; kwargs...)
    @error("score not implemented for sufficient statistics type $(typeof(ss)) \n
        you may need to implement a custom sufficient statistics type")
end

### ========================================================================================

struct CategoricalSuffStats{M, T} <: SuffStats
    h::SVector{M, T}
    n::Int
end

function CategoricalSuffStats(num_categories::Int)
    h = SVector{num_categories, Int}(zeros(Int, num_categories))
    return CategoricalSuffStats{num_categories, Int}(h, 0)
end

@inline function add_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] += 1
    ss = @set ss.n += 1
    return ss
end

@inline function add_sample(ss::CategoricalSuffStats, ::Nothing)
    @reset ss.n += 1
    return ss
end

@inline function remove_sample(ss::CategoricalSuffStats, sample::Int)
    ss = @set ss.h[sample] -= 1
    ss = @set ss.n -= 1
    return ss
end

@inline function remove_sample(ss::CategoricalSuffStats, ::Nothing)
    @reset ss.n -= 1
    return ss
end

function make_k_block(k, ::Val{:categorical}; num_categories, kwargs...)
    k_block = SymArray{CategoricalSuffStats{num_categories, Int}}(undef, k, k)
    fill!(k_block, CategoricalSuffStats(num_categories))
    return k_block
end

@inline function score(ss::CategoricalSuffStats; kwargs...)
    n = max(ss.n, 1)
    return n - sum(abs2, ss.h) / n
end

### ========================================================================================

struct BernoulliSuffStats{T} <: SuffStats
    h::T
    n::T
end

function BernoulliSuffStats()
    return BernoulliSuffStats{Int}(0, 0)
end

function add_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h += 1)
    @reset ss.n += 1
    return ss
end

function add_sample(ss::BernoulliSuffStats, ::Nothing)
    @reset ss.n += 1
    return ss
end

function remove_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h -= 1)
    @reset ss.n -= 1
    return ss
end

function remove_sample(ss::BernoulliSuffStats, ::Nothing)
    @reset ss.n -= 1
    return ss
end

function make_k_block(k, ::Val{:binary}; kwargs...)
    k_block = SymArray{BernoulliSuffStats{Int}}(undef, k, k)
    fill!(k_block, BernoulliSuffStats())
    return k_block
end

@inline function score(ss::BernoulliSuffStats; kwargs...)
    n = max(ss.n, 1)
    p = ss.h / n
    return n * (xlogx(1 - p) + xlogx(p))
end

### ========================================================================================

struct GenericSuffStats{T} <: SuffStats
    samples::Vector{T}
end

function GenericSuffStats{T}() where {T}
    return GenericSuffStats{T}(Vector{T}())
end

function add_sample(ss::GenericSuffStats, sample)
    append!(ss.samples, sample)
    return ss
end

function remove_sample(ss::GenericSuffStats, sample)
    index = findfirst(==(sample), ss.samples)
    if index !== nothing
        deleteat!(ss.samples, index)
    end
    return ss
end

function make_k_block(k, generic; data::AbstractArray, kwargs...)
    @warn "Using GenericSuffStats may lead to high memory usage for large datasets.
         Consider using more specialized sufficient statistics types when possible."
    k_block = SymArray{GenericSuffStats{eltype(data)}}(undef, k, k)
    for j in 1:k, i in 1:k
        k_block[i, j] = GenericSuffStats{eltype(data)}()
    end
    return k_block
end

function score(ss::GenericSuffStats; dist::D, kwargs...) where {D}
    if dist === nothing
        @error("No distribution provided for scoring GenericSuffStats")
    end
    d = fit(D, ss.samples)
    return -sum(logpdf.(Ref(d), ss.samples))
end

### ========================================================================================

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
            es.block_ss[gi, gj] = add_sample(es.block_ss[gi, gj], edge_value)
            es.block_ss_swap[gi, gj] = add_sample(es.block_ss_swap[gi, gj], edge_value)
        end
    end
end

# TODO: allow for non-symmetric data
@inline function score(matrix_ss::SymArray, data, node_labels; dist = nothing, norm = 1.0)
    total_loss = 0.0
    for m in matrix_ss.uppertrian.nzval
        total_loss += score(m; dist = dist, data = data, node_labels = node_labels)
    end
    # @inbounds for j in axes(matrix_ss, 2)
    #     for i in 1:j
    #         inter = score(
    #             matrix_ss[i, j]; dist = dist, data = data, node_labels = node_labels)
    #         total_loss += inter
    #     end
    # end
    return total_loss / norm
end

function GreedySuffStats(
        data, node_labels; type_suff_stats = :categorical, max_iter = 10000,
        node_swap_rule = RandomGroupSwap(), stop_rule = PreviousBestValue(5_000, Inf, :min),
        kwargs...)
    # derive user input
    k = length(unique(node_labels))

    # allocate sufficient statistics blocks
    block_ss = make_k_block(k, Val(type_suff_stats); data = data, kwargs...)
    block_ss_swap = make_k_block(k, Val(type_suff_stats); data = data, kwargs...)

    # create estimator
    return GreedySuffStats{typeof(block_ss), typeof(node_swap_rule), typeof(stop_rule)}(
        block_ss, block_ss_swap, node_swap_rule, stop_rule, max_iter)
    return es
end

function estimate(
        es::GreedySuffStats,
        data,
        node_labels_init;
        progress = true,
        dist = nothing,
        iter_progress = 5000
)
    # Initialize node labels
    node_labels = copy(node_labels_init)
    n = length(node_labels)
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

    current_loss = score(es.block_ss, data, node_labels, dist = dist, norm = n_edges)
    es.stop_rule.previous_best_value = current_loss
    # Main optimization loop
    for iter in 1:(es.max_iter)
        # Select two nodes to potentially swap
        index1, index2 = select_indices_swap(node_labels, es.node_swap_rule)

        group1 = node_labels[index1]
        group2 = node_labels[index2]

        @inbounds for j in axes(data, 2)

            # this check is slow ! (+ 6 μs per iteration on n=2000)
            if j == index1 || j == index2
                continue
            end
            # extract data
            groupj = node_labels[j]
            edge_value_1 = data[j, index1]
            edge_value_2 = data[j, index2]

            es.block_ss_swap[group1, groupj] = remove_sample(
                es.block_ss_swap[group1, groupj], edge_value_1)
            es.block_ss_swap[group2, groupj] = add_sample(
                es.block_ss_swap[group2, groupj], edge_value_1)

            es.block_ss_swap[group2, groupj] = remove_sample(
                es.block_ss_swap[group2, groupj], edge_value_2)
            es.block_ss_swap[group1, groupj] = add_sample(
                es.block_ss_swap[group1, groupj], edge_value_2)
        end

        # tentative swap
        @inbounds node_labels[index1], node_labels[index2] = group2, group1
        new_loss = score(es.block_ss_swap, data, node_labels, dist = dist, norm = n_edges)

        if new_loss < current_loss
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

    return node_labels
end
