struct Assignment{T <: Real}
    number_nodes::Int
    number_groups::Int
    group_size::Int
    proportion::T

    node_labels::Vector{Int}
    counts::Matrix{Int}
    realized::Matrix{Float64}
    estimated_theta::Matrix{Float64}

    likelihood::Float64

    function Assignment(A, node_labels, h)
        number_groups = length(unique(node_labels))
        number_nodes = length(node_labels)
        group_size = number_nodes % h
        estimated_theta = zeros(Float64, number_nodes, number_groups)
        counts = zeros(Int64, number_groups, number_groups)
        realized = zeros(Int64, number_groups, number_groups)
        size_groups = [count(==(element), node_labels) for element in unique(node_labels)]

        @inbounds @simd for k in 1:number_groups
            for l in k:number_groups
                realized[k, l] = sum(A[node_labels .== k, node_labels .== l])
                realized[l, k] = realized[k, l]
                counts[k, l] = size_groups[k] * size_groups[l]
                counts[l, k] = counts[k, l]
            end
        end

        @inbounds @simd for k in 1:number_groups
            counts[k, k] = size_groups[k] * (size_groups[k] - 1) / 2
            realized[k, k] = sum(A[node_labels .== k, node_labels .== k]) / 2
        end

        estimated_theta = realized ./ counts
        likelihood = compute_log_likelihood(number_groups, estimated_theta, counts,
                                            number_nodes)

        new(number_nodes,
            number_groups,
            node_labels,
            h,
            group_size,
            counts,
            realized,
            estimated_theta,
            likelihood)
    end

    function Assignment(a::Assignment, likelihood)
        new(a.number_nodes,
            a.number_groups,
            a.node_labels,
            a.group_size,
            a.counts,
            a.realized,
            a.estimated_theta,
            likelihood)
    end
end

function initialize(A, h::Int, optimizer)
    node_labels = initialise_node_labels(A, h)
    old_store = Assignment(A, node_labels, h)
    new_store = deepcopy(oldstore)
    history = MVHistory([
                            :likelihood => QHistory(Float64),
                            :best_likelihood => QHistory(Float64),
                            :proposal_likelihood => QHistory(Float64),
                        ])
    return old_store, new_store, history
end

function initialise_node_labels(A, h::Int)
    error("Not yet implemented")
end

function compute_log_likelihood(number_groups, estimated_theta, counts, number_nodes)
    loglik = 0.0
    @inbounds @simd for i in 1:number_groups
        for j in i:number_groups
            θ = estimated_theta[i, j]
            θ_c = θ < 0 ? 1e-16 : (θ > 1 ? 1 - 1e-16 : θ)
            loglik += (θ_c * log(θ_c) + (1 - θ_c) * log(1 - θ_c)) * counts[i, j]
        end
    end
    return loglik / number_nodes
end



function deepcopy!(a::Assignment, b::Assignment)
    a.number_nodes = b.number_nodes
    a.number_groups = b.number_groups
    a.node_labels .= b.node_labels
    a.group_size = b.group_size
    a.counts .= b.counts
    a.realized .= b.realized
    a.estimated_theta .= b.estimated_theta
    a.likelihood = b.likelihood
end
