mutable struct Assignment{T}
    const group_size::GroupSize{T}

    const node_labels::Vector{Int}
    const counts::Matrix{Int}
    const realized::Matrix{Float64}
    const estimated_theta::Matrix{Float64}

    likelihood::Float64

    function Assignment(A, node_labels, group_size::GroupSize{T}) where {T}
        number_groups = length(group_size)
        estimated_theta = zeros(Float64, size(A, 1), number_groups)

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
                                            size(A, 1))

        new{T}(group_size,
               node_labels,
               counts,
               realized,
               estimated_theta,
               likelihood)
    end
end

"""
    compute_log_likelihood(number_groups, estimated_theta, counts, number_nodes)

Compute the scaled log-likelihood in terms of communities:
```math
l(z;A) = \\frac{1}{n} \\sum_{g_1 = 1}^{G} \\sum_{g_2 \\geq g_1}^{G}  \\left[ \\theta_{g_1g_2} \\log(\\theta_{g_1g_2}) + (1 - \\theta_{g_1g_2}) \\log(1 - \\theta_{g_1g_2}) \\right] \\cdot c_{g_1g_2},
```

where ``c_{g_1g_2}`` (``\\theta_{g_1g_2}``) is the number of possible edges (estimated
probability) between communities ``g_1`` and ``g_2``, ``n`` is the number of nodes, and
``z_i ∈ \\{1, \\dots, G\\}`` is the community assignment of node ``i``.
"""
function compute_log_likelihood(number_groups, estimated_theta, counts, number_nodes)
    loglik = 0.0
    @inbounds @simd for i in 1:number_groups
        for j in i:number_groups
            θ = estimated_theta[i, j]
            θ_c = θ <= 0 ? 1e-14 : (θ >= 1 ? 1 - 1e-14 : θ)
            loglik += (θ_c * log(θ_c) + (1 - θ_c) * log(1 - θ_c)) * counts[i, j]
        end
    end
    return loglik / number_nodes
end

"""
    compute_log_likelihood(assignment::Assignment)

Compute the scaled log-likelihood of the assignment.

```math
	l(z;A) = \\frac{1}{n}\\sum\\limits_{i=1}^n \\sum\\limits_{j>i}^n  \\left[ A_{ij} \\log(\\hat{\\theta}_{z_i z_j}) + (1 - A_{ij}) \\log(1 - \\hat{\\theta}_{z_i z_j}) \\right],
```

where ``\\hat{\\theta}_{ab}`` is the estimated probability of an edge between communities
``a`` and ``b``

```math
    \\hat{\\theta}_{ab} = \\frac{\\sum\\limits_{i<j} A_{ij} \\mathbb{1}(z_i = a, z_j = b) }{\\sum\\limits_{i<j} \\mathbb{1}(z_i = a, z_j = b)}.
```
"""
function compute_log_likelihood(assignment::Assignment)
    compute_log_likelihood(length(assignment.group_size),
                           assignment.estimated_theta,
                           assignment.counts,
                           sum(assignment.group_size))
end

function deepcopy!(a::Assignment, b::Assignment)
    a.node_labels .= b.node_labels
    a.counts .= b.counts
    a.realized .= b.realized
    a.estimated_theta .= b.estimated_theta
    a.likelihood = b.likelihood
end
