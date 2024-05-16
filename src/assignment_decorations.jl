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


"""
    deepcopy!(a::Assignment, b::Assignment)

Deep copy the attributes of `b` to `a`.
"""
function deepcopy!(a::Assignments, b::Assignments)
    a.node_labels .= b.node_labels
    a.estimated_theta .= b.estimated_theta
    a.likelihood = b.likelihood
    a.edge_labels .= b.edge_labels
end


get_parameters(d::Union{UnivariateDistribution, MultivariateDistribution}) = params(d)
fit(d::Union{UnivariateDistribution, MultivariateDistribution}, x) = fit_mle(typeof(d), x)
fit(d::Categorical, x) = fit_mle(typeof(d), length(d.p), x)



function update(a::Assignments{T,D}, A::AbstractArray{S}, s) where {T,D,S}
    new_node_labels = copy(a.node_labels)
    if a.node_labels[s[1]] != a.node_labels[s[2]]
        new_node_labels[s[1]] = a.node_labels[s[2]]
        new_node_labels[s[2]] = a.node_labels[s[1]]
        return Assignments(A, new_node_labels, a.estimated_theta[1])
    else
        return a
    end
end


function initialize(A,d, h, starting_assignment_rule, record_trace)
    node_labels, group_size = initialize_node_labels(A, h, starting_assignment_rule)
    proposal = Assignments(A, node_labels, d, group_size)
    current = deepcopy(proposal)
    best = deepcopy(proposal)
    history = initialize_history(best, current, proposal, record_trace)
    return best, current, proposal, history, A
end


function map_to_old_format(a::Assignments)
    return map_to_old_format(a.estimated_theta, a.map, length(a.group_size))
end

function map_to_old_format(dists, map, n_groups)
    t = zeros(n_groups, n_groups)
    for i in 1:n_groups
        for j in 1:n_groups
            tuple_ordered = i ≤ j ? (i, j) : (j, i)
            t[i, j] = params(dists[map[tuple_ordered]])[1]
        end
    end
    return t
end


function select_bandwidth(A, d)::Int
    h = oracle_bandwidth(A,d)
    return max(2, min(size(A)[1], round(Int, h)))
end

# idea is to take the least regular coordinate of the decorated graphon
# and use it for setting the bandwidth → we take the minimum α in a sense
function oracle_bandwidth(A::Matrix{S}, d) where {S}
    K = unique(A)
    h = size(A,1)÷2
    A_inter = BitMatrix(undef, size(A)...)
    for k in K
        indices = findall( x -> x == k, A)
        A_inter .= 0
        A_inter[indices] .= 1
        h = min(h, oracle_bandwidth(A_inter))
    end
    return h
end


function from_adj_to_lower_triag_with_diag(A::AbstractMatrix)
    return A[CartesianIndex.((i, j) for j in axes(A, 2) for i in axes(A, 1) if i ≥ j)]
end

function from_lower_triag_with_diag_to_adj(B)
    n = Int(sqrt(2 * length(B) + 0.25) - 0.5)
    A = similar(B, n, n)
    A[CartesianIndex.((i, j) for j in axes(A, 2) for i in axes(A, 1) if i ≥ j)] .= B
    A = A + A'
    for i in 1:n
        A[i, i] = A[i, i] .÷ 2
    end
    return A
end
