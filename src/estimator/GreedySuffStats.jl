abstract type SuffStats end

function add_sample end
function remove_sample end

add_sample(suffstats::SuffStats, sample, i, j) = add_sample(suffstats, sample)
remove_sample(suffstats::SuffStats, sample, i, j) = remove_sample(suffstats, sample)
function make_k_block end
function score end

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
    return ss.n - sum(abs2, ss.h) / max(ss.n, 1)
end

### ========================================================================================

struct BernoulliSuffStats{T} <: SuffStats
    h::T
    n::T
end

function BernoulliSuffStats()
    return BernoulliSuffStats{Int}(0, 0)
end

@inline function add_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h += 1)
    @reset ss.n += 1
    return ss
end

@inline function add_sample(ss::BernoulliSuffStats, ::Nothing)
    @reset ss.n += 1
    return ss
end

@inline function remove_sample(ss::BernoulliSuffStats, sample::Bool)
    sample && (@reset ss.h -= 1)
    @reset ss.n -= 1
    return ss
end

@inline function remove_sample(ss::BernoulliSuffStats, ::Nothing)
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

abstract type GenericSuffStatsType <: SuffStats end

struct GenericSuffStats{T} <: GenericSuffStatsType
    samples::Vector{T}
end

function GenericSuffStats(::AbstractArray{T}) where {T}
    return GenericSuffStats{T}(Vector{T}())
end

function get_samples(ss::GenericSuffStats)
    return ss.samples
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
        k_block[i, j] = GenericSuffStats(data)
    end
    return k_block
end

# use indices rather than pushing and deleting samples for better performance ?
# struct GenericSuffStatsIndex{T} <: GenericSuffStatsType
#     indices::Vector{Tuple{Int, Int}}
#     data::T
# end

# function get_samples(ss::GenericSuffStatsIndex)
#     return [ss.data[i, j] for (i, j) in ss.indices]
# end

# function GenericSuffStatsIndex{T}(data::T) where {T}
#     return GenericSuffStatsIndex{T}(Vector{Tuple{Int, Int}}(), data)
# end

# function add_sample(ss::GenericSuffStatsIndex, sample, i, j)
#     push!(ss.indices, (i, j))
#     return ss
# end

function score(ss::GenericSuffStatsType; dist::D, kwargs...) where {D}
    if dist === nothing
        @error("No distribution provided for scoring GenericSuffStats")
    end
    samples = get_samples(ss)
    d = fit(D, samples)
    return -sum(logpdf.(d, samples))
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
            es.block_ss[gi, gj] = add_sample(es.block_ss[gi, gj], edge_value, i, j)
            es.block_ss_swap[gi, gj] = add_sample(
                es.block_ss_swap[gi, gj], edge_value, i, j)
        end
    end
end

# TODO: allow for non-symmetric data
@inline function score(matrix_ss::SymArray, data, node_labels; dist = nothing, norm = 1.0)
    total_loss = 0.0
    for m in matrix_ss.uppertrian.nzval
        total_loss += score(m; dist = dist, data = data, node_labels = node_labels)
    end
    return total_loss / norm
end

@inline function score(matrix_ss, data, node_labels; dist = nothing, norm = 1.0)
    total_loss = 0.0
    @inbounds for j in axes(matrix_ss, 2)
        for i in 1:j
            inter = score(
                matrix_ss[i, j]; dist = dist, data = data, node_labels = node_labels)
            total_loss += inter
        end
    end
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

function estimate!(
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
