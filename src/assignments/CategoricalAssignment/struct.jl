mutable struct CategoricalData{F}
    counts::Matrix{Int}
    realized::Matrix{Array{Int,2}}
    estimated_theta::Matrix{Array{F,2}}
    A::Matrix{Int}
    log_likelihood::F
end

const CategoricalAssignment{T, F} = Assignment{T, CategoricalData{F}}
const CategoricalInitRule{S, F} = InitRule{S, Val{CategoricalData}}

function CategoricalAssignment(
        g, group_size::GroupSize, node_labels::Vector{Int})
    categorical_data = make_categorical_data(g, node_labels, group_size)
    return Assignment(group_size, node_labels, categorical_data)
end

function make_assignment(g, h, init_rule::CategoricalInitRule)
    group_size,
    node_labels = initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)
    return CategoricalAssignment(g, group_size, node_labels)
end

function make_categorical_data(g, node_labels, group_size)
    number_groups = length(group_size)
    n = length(node_labels)
    A, num_categories = categorical_matrix(g)
    counts = zeros(Int, number_groups, number_groups)
    realized = zeros(Int, number_groups, number_groups, num_categories)
    @inbounds @simd for k in 1:number_groups
        for l in k:number_groups
            for m in 1:num_categories
                if k == l
                    c = group_size[k] * (group_size[k] - 1) ÷ 2
                    r = sum(A[node_labels .== k, node_labels .== l] .== m)/2
                else
                    c = group_size[k] * group_size[l]
                    r = sum(A[node_labels .== k, node_labels .== l] .== m)
                end
                realized[k, l, m] = r
                realized[l, k, m] = realized[k, l, m]
                counts[k, l] = c
                counts[l, k] = c
            end
        end
    end

    estimated_theta = realized ./ counts
    ll = compute_log_likelihood(estimated_theta, counts)
    return CategoricalData(counts, realized, estimated_theta, A, ll)
end



function compute_log_likelihood_1(estimated_theta, counts)
    loglik = 0.0
    @inbounds @simd for coord in CartesianIndices(estimated_theta)
        if coord[2] <= coord[1]
            θ = estimated_theta[coord]
            θ_c = θ <= 0 ? 1e-14 : (θ >= 1 ? 1 - 1e-14 : θ)
            loglik += (θ_c * log(θ_c) + (1 - θ_c) * log(1 - θ_c)) *
                      counts[coord]
        end
    end
    return loglik
end
