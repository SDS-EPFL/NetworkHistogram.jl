struct Assignment{T}
    group_size::GroupSize{T}

    node_labels::Vector{Int}
    counts::Matrix{Int}
    realized::Matrix{Float64}
    estimated_theta::Matrix{Float64}

    likelihood::Float64

    function Assignment(A, node_labels, group_size::GroupSize{T}) where {T}
        number_groups = length(group_size)
        estimated_theta = zeros(Float64, number_nodes, number_groups)

        counts = zeros(Int64, number_groups, number_groups)
        realized = zeros(Int64, number_groups, number_groups)

        @inbounds @simd for k in 1:number_groups
            for l in k:number_groups
                realized[k, l] = sum(A[node_labels .== k, node_labels .== l])
                realized[l, k] = realized[k, l]
                counts[k, l] = group_size[k] * group_size[l]
                counts[l, k] = counts[k, l]
            end
        end

        @inbounds @simd for k in 1:number_groups
            counts[k, k] = group_size[k] * (group_size[k] - 1) ÷ 2
            realized[k, k] = sum(A[node_labels .== k, node_labels .== k]) ÷ 2
        end

        estimated_theta = realized ./ counts
        likelihood = compute_log_likelihood(number_groups, estimated_theta, counts,
                                            number_nodes)

        new{T}(
            group_size,
            node_labels,
            counts,
            realized,
            estimated_theta,
            likelihood
        )
    end

    function Assignment(a::Assignment{T}, likelihood) where {T}
        new{T}(
            a.group_size,
            a.node_labels,
            a.counts,
            a.realized,
            a.estimated_theta,
            likelihood
        )
    end
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
    a.node_labels .= b.node_labels
    a.counts .= b.counts
    a.realized .= b.realized
    a.estimated_theta .= b.estimated_theta
    return Assignment(a, b.likelihood)
end
