mutable struct Assignments{T,D}
    const group_size::GroupSize{T}
    const node_labels::Vector{Int}
    const edge_labels::Vector{Int}
    const estimated_theta::Array{D}
    const map::Dict{Tuple{Int, Int}, Int}
    likelihood::Float64

    function Assignments(A::AbstractArray{S,2}, node_labels, d::D, group_size::GroupSize{T} = GroupSize(node_labels)) where {D, S, T}
        edge_indices = lowertriag(A)
        edge_labels, map = node_labels_to_edge_labels(node_labels, edge_indices)
        estimated_theta, ll = estimate_theta(A, edge_labels, d, edge_indices)
        new{T, D}(
            group_size, node_labels, edge_labels, estimated_theta, map,ll)
    end
end



#### TODO: implement with weights → no assighment of edge variables, fit weighted
# mixture ?


function estimate_theta(A::AbstractArray{S}, edge_labels, d::D, edge_indices) where {S,D}
    k_edge_groups = length(unique(edge_labels))
    estimated_theta = Array{D}(undef, k_edge_groups)
    ll = 0.0
    for edge_group in 1:k_edge_groups
        X = A[edge_indices[edge_labels .== edge_group]]
        estimated_theta[edge_group] = fit(d, X)
        ll += compute_log_likelihood(X, estimated_theta[edge_group])
    end
    return estimated_theta,ll
end



function compute_log_likelihood(
        X::AbstractVector{S}, estimated_theta::D) where {S, D}
    ll = 0.0
    for edge in X
        ll += logdensityof(estimated_theta, edge)
    end
    return ll
end



function compute_log_likelihood(X::AbstractVector{S}, edge_labels, estimated_theta::Vector{D}) where {S,D}
    ll = 0.0
    for edge_index in eachindex(X)
        edge_group = edge_labels[edge_index]
        ll += logdensityof(estimated_theta[edge_group], X[edge_index])
    end
    return ll
   i
end


function lowertriag(A::AbstractMatrix)
    CartesianIndex.((i, j) for j in axes(A, 2) for i in axes(A, 1) if i > j)
end


function node_labels_to_edge_labels(node_labels::Vector{Int}, edge_indices)
    edge_labels = zeros(Int, length(edge_indices))
    map = Dict{Tuple{Int,Int},Int}()
    k = 1
    for i in 1:length(unique(node_labels))
        for j in i:length(unique(node_labels))
            map[(i, j)] = k
            k += 1
        end
    end
    for (index, c) in enumerate(edge_indices)
        i, j = Tuple(c)
        node_tuple = node_labels[i] ≤ node_labels[j] ? (node_labels[i], node_labels[j]) : (node_labels[j], node_labels[i])
        edge_labels[index] = node_tuple |> x -> map[x]
    end
    return edge_labels, map
end



get_parameters(d::Union{UnivariateDistribution, MultivariateDistribution}) = params(d)
fit(d::Union{UnivariateDistribution, MultivariateDistribution}, x) = fit_mle(typeof(d), x)
fit(d::Categorical, x) = fit_mle(typeof(d), length(d.p), x)
