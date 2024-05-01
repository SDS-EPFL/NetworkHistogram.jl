mutable struct Assignment{T,D}
    const group_size::GroupSize{T}
    const node_labels::Vector{Int}
    const edge_labels::Vector{Int}
    const estimated_theta::Array{D}
    likelihood::Float64

    function Assignment(A::Matrix{S}, node_labels, d::D) where {D, S}
        group_size = GroupSize(node_labels)
        edge_labels, map = node_labels_to_edge_labels(node_labels)
        println(edge_labels)
        X = adjacency_matrix_to_ordered_edge_list(A)
        estimated_theta, ll = estimate_theta(X, edge_labels, d)
        new{typeof(group_size[1]), D}(
            group_size, node_labels, edge_labels, estimated_theta, ll)
    end
end






function estimate_theta(X::Array{S}, edge_labels, d::D) where {S,D<:UnivariateDistribution}
    k_edge_groups = length(unique(edge_labels))
    estimated_theta = Array{D}(undef, k_edge_groups)
    for edge_group in eachindex(estimated_theta)
        edges = findall(edge_labels .== edge_group)
        estimated_theta[edge_group] = fit(d, X[edges])
    end
    ll = compute_log_likelihood(X, edge_labels, estimated_theta)
    return estimated_theta,ll
end


function compute_log_likelihood(X, edge_labels, estimated_theta)
    ll = 0.0
    for edge_index in eachindex(X)
        edge_group = edge_labels[edge_index]
        ll += logdensityof(estimated_theta[edge_group], X[edge_index])
    end
    return ll

end


function adjacency_matrix_to_ordered_edge_list(A::Matrix{S}) where {S}
    n = size(A,1)
    edge_list = Array{S,1}(undef, n*(n-1)÷2)
    index = 1
    for i in 1:n
        for j in i+1:n
            edge_list[index] = A[i,j]
            index += 1
        end
    end
    return edge_list
end




function node_labels_to_edge_labels(node_labels::Vector{Int})
    n = length(node_labels)
    k = maximum(node_labels)
    map = Dict{Tuple{Int,Int},Int}()
    edge_labels = zeros(Int, n*(n-1)÷2)
    index = 1
    index_group = 1
    for i in 1:n
        for j in i+1:n
            if !haskey(map, (node_labels[i], node_labels[j]))
                map[(node_labels[i], node_labels[j])] = index_group
                index_group += 1
            end
            edge_labels[index] = map[(node_labels[i], node_labels[j])]
            index += 1
        end
    end
    return edge_labels, map
end



get_parameters(d::Union{UnivariateDistribution, MultivariateDistribution}) = params(d)
fit(d::Union{UnivariateDistribution, MultivariateDistribution}, x) = fit_mle(typeof(d), x)
fit(d::Categorical, x) = fit_mle(typeof(d), length(d.p), x)
