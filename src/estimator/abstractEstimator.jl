abstract type SBMEstimator end

struct SumGreedyEstimator{C, S, P} <: SBMEstimator
    counts::C # counts of possible edges between groups
    counts_swap::C
    realized::S # sums of observed Λ between groups
    realized_swap::S
    max_iter::Int
    stop_rule::P
end

function init!(es::SumGreedyEstimator, data, initial_labels)
    for j in axes(data, 2)
        label_j = initial_labels[j]
        for i in axes(data, 1) # double counting edges in undirected graphs
            if !isnothing(data[i, j]) && i < j
                add_counts!(es.realized[initial_labels[i], label_j], data[i, j])
                add_counts!(es.realized_swap[initial_labels[i], label_j], data[i, j])
                es.counts[initial_labels[i], label_j] += 1
                es.counts_swap[initial_labels[i], label_j] += 1
            end
        end
    end
end

function estimate(es::SumGreedyEstimator, data, initial_labels)
    init!(es, data, initial_labels)
    loss = loss_function(es.realized, es.counts)
    es.stop_rule.previous_best_value = loss
    ## optim
    node_labels = copy(initial_labels)
    k = length(unique(node_labels))
    pbar = ProgressUnknown(enabled = true, showspeed = true, desc = "Greedy search: ")
    for iter in 1:(es.max_iter)
        next!(pbar)
        groups = StatsBase.sample(1:k, 2; replace = false)
        index1 = rand(findall(x -> x == groups[1], node_labels))
        index2 = rand(findall(x -> x == groups[2], node_labels))
        #index1, index2 = StatsBase.sample(1:length(node_labels), 2; replace = false)
        g1 = node_labels[index1]
        g2 = node_labels[index2]
        if g1 == g2
            continue
        end
        edges_index1 = view(data, :, index1)
        edges_index2 = view(data, :, index2)
        for j in axes(data, 1)
            if j == index1 || j == index2
                continue
            end
            gj = node_labels[j]
            remove_counts!(es.realized_swap[g1, gj], edges_index1[j])
            es.counts_swap[g1, gj] -= 1
            add_counts!(es.realized_swap[g2, gj], edges_index1[j])
            es.counts_swap[g2, gj] += 1
            remove_counts!(es.realized_swap[g2, gj], edges_index2[j])
            es.counts_swap[g2, gj] -= 1
            add_counts!(es.realized_swap[g1, gj], edges_index2[j])
            es.counts_swap[g1, gj] += 1
        end
        node_labels[index1] = g2
        node_labels[index2] = g1

        loss_new = loss_function(es.realized_swap, es.counts_swap)
        if loss_new < es.stop_rule.previous_best_value
            es.stop_rule.previous_best_value = loss_new
            es.stop_rule.iterations_since_best = 0
            deepcopy!(es.realized, es.realized_swap)
            copy!(es.counts, es.counts_swap)
            loss = loss_new
        else
            # revert swap
            node_labels[index1] = g1
            node_labels[index2] = g2
            deepcopy!(es.realized_swap, es.realized)
            copy!(es.counts_swap, es.counts)
            es.stop_rule.iterations_since_best += 1
        end

        if es.stop_rule.iterations_since_best >= es.stop_rule.k
            @info "Stopping criterion met at iteration $iter"
            finish!(pbar)
            break
        end
    end
    return node_labels, losses
end

function loss_function(realized, counts)
    loss = 0.0
    @inbounds for j in axes(realized, 2)
        for i in axes(realized, 1)
            if i <= j
                # θ = realized[i, j] ./ counts[i, j]
                # loss += sum(xlogx, θ) * counts[i, j]
                loss += sum(counts[i, j] - sum(abs2, realized[i, j]) / counts[i, j])
            end
        end
    end
    return loss
end

function add_counts!(parameter::AbstractArray, data_value::AbstractArray)
    @inbounds parameter .+= data_value
end

function remove_counts!(parameter::AbstractArray, data_value::AbstractArray)
    @inbounds parameter .-= data_value
end

function add_counts!(parameter::AbstractArray, data_value::Real)
    @inbounds parameter[data_value] += 1
end

function remove_counts!(parameter::AbstractArray, data_value::Real)
    @inbounds parameter[data_value] -= 1
end

##
using Makie

function Makie.convert_arguments(
        ::Type{<:AbstractPlot}, graphon::DecoratedSBM, k::Int = 1)
    x = collect(0:0.01:1)
    return (x, x, [_extract_param(graphon(xi, yi), k) for xi in x, yi in x])
end

function _extract_param(d::Distribution, k::Int)
    return params(d)[k]
end

function _extract_param(d::DiscreteNonParametric, k::Int)
    return params(d)[2][k]
end
