"""
    Assignment{E, D, F, W, V <: AbstractVector{Int}}

Represents a network histogram: a partition of nodes into groups along with
edge distributions between groups.

# Fields
- `node_labels::V`: Vector assigning each node to a group (1-indexed)
- `edges::EdgeList{E}`: The observed edge data
- `dists::EdgeList{D}`: Fitted distributions for each edge
- `θ::SymArray{D}`: Symmetric matrix of aggregated distributions between groups
- `log_likelihood::SymArray{F}`: Symmetric matrix of log-likelihoods for each group pair
- `additional_workspace::W`: Optional workspace for optimization algorithms

# Type Parameters
- `E`: Type of edge observations
- `D`: Type of fitted distributions
- `F`: Type for log-likelihood values (typically Float64)
- `W`: Type for additional workspace data
- `V`: Vector type for node labels

# Examples
```julia
# Create assignment from node labels and edge data
node_labels = [1, 1, 2, 2, 3]
edges = EdgeList(adjacency_matrix)
dist = Dist(Bernoulli(0.5))
assignment = Assignment(node_labels, edges, dist)

# Query assignment properties
k = number_groups(assignment)
n = number_nodes(assignment)
ll = loglikelihood(assignment)
group_i = group(assignment, node_i)
```

See also: [`BlockModel`](@ref), [`EdgeList`](@ref), [`Dist`](@ref)
"""
mutable struct Assignment{E, D, F, W, V <: AbstractVector{Int}}
    node_labels::V
    const edges::EdgeList{E}
    const dists::EdgeList{D}
    θ::SymArray{D}
    log_likelihood::SymArray{F}
    additional_workspace::W
end

"""
    number_nodes(a::Assignment)

Return the number of nodes in the network.
"""
@inline number_nodes(a::Assignment) = length(a.node_labels)

"""
    number_groups(a::Assignment)

Return the number of groups (blocks) in the partition.
"""
@inline number_groups(a::Assignment) = size(a.θ, 1)

"""
    proportions(a::Assignment)

Calculate the proportion of nodes in each group.

# Returns
- Vector of proportions summing to 1.0
"""
function proportions(a::Assignment)
    return counts(a.node_labels) / number_nodes(a)
end

"""
    loglikelihood(a::Assignment)

Calculate the total log-likelihood of the assignment.

The log-likelihood measures how well the stochastic block model (with the current
node partition) fits the observed network data.

# Returns
- `Float64`: Total log-likelihood value
"""
@inline function loglikelihood(a::Assignment)
    return FastSymArray.sum_tri_with_diag(a.log_likelihood)
end

"""
    group(a::Assignment, node::Int)

Get the group label for a specific node.

# Arguments
- `a::Assignment`: The assignment
- `node::Int`: Node index (1-indexed)

# Returns
- `Int`: Group index that the node belongs to
"""
@inline function group(a::Assignment, node::Int)
    @boundscheck checkbounds(a.node_labels, node)
    @inbounds return a.node_labels[node]
end

"""
    get_edges_in_groups(a::Assignment, g1::Int, g2::Int)

Extract all edges between two groups.

# Arguments
- `a::Assignment`: The assignment
- `g1::Int`: First group index
- `g2::Int`: Second group index

# Returns
- `Vector{E}`: Vector of edge values between the two groups

# Note
For within-group edges (g1 == g2), only returns edges where i < j to avoid duplicates.
"""
function get_edges_in_groups(a::Assignment, g1::Int, g2::Int)
    return get_edges_in_groups(a.node_labels, a.edges, g1, g2)
end

function get_edges_in_groups(node_labels, edges_all, g1, g2)
    edges = Vector{edge_type(edges_all)}()
    nodes_g1 = findall(x -> x == g1, node_labels)
    nodes_g2 = findall(x -> x == g2, node_labels)

    for u in nodes_g1
        for (v, e) in iterate_neighbors(edges_all, u)
            if v in nodes_g2 && ((g1 == g2 && u < v) || g1 != g2)
                push!(edges, e)
            end
        end
    end
    return edges
end

"""
    Assignment(node_labels, edge_list::EdgeList{E}, dist::Dist{D}) where {E, D}

Construct an Assignment from node labels, edge data, and a reference distribution.

This constructor fits the distribution to the data, computes the block-level parameters
θ, and calculates the log-likelihood.

# Arguments
- `node_labels`: Vector of group assignments for each node
- `edge_list::EdgeList{E}`: Edge observations
- `dist::Dist{D}`: Reference distribution to fit to the data

# Example
```julia
node_labels = [1, 1, 2, 2]
edges = EdgeList(A)
dist = Dist(Bernoulli(0.5))
assignment = Assignment(node_labels, edges, dist)
```
"""
function Assignment(
        node_labels, edge_list::EdgeList{E}, dist::Dist{D}) where {E, D}
    dists = fit(dist, edge_list)
    θ, ll = _compute_theta_and_ll(node_labels, dists, edge_list, dist)
    return Assignment(node_labels, edge_list, dists, θ, ll, nothing)
end

# Internal function to compute θ parameters and log-likelihood for each group pair
function _compute_theta_and_ll(node_labels, dists::EdgeList{Dist{D}},
        edge_list::EdgeList{E}, dist::Dist{D}) where {E, D}
    number_groups = length(unique(node_labels))
    θ = SymArray(number_groups, zero(dist))
    log_likelihood = SymArray(number_groups, 0.0)

    # Aggregate distributions for each group pair
    for u in 1:nodes(dists)
        g1 = node_labels[u]
        for (v, d) in iterate_neighbors(dists, u)
            g2 = node_labels[v]
            if u < v
                θ[g1, g2] = add_to(θ[g1, g2], d)
            end
        end
    end

    # Compute log-likelihood for each group pair
    for u in 1:nodes(dists)
        g1 = node_labels[u]
        for (v, e) in iterate_neighbors(edge_list, u)
            g2 = node_labels[v]
            if u > v
                log_likelihood[g1, g2] += logpdf(θ[g1, g2], e)
            else
                break
            end
        end
    end
    return θ, log_likelihood
end
