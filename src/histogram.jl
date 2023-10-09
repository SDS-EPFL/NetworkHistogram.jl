struct GraphHist{T,M}
    θ::Array{T,M}
    node_labels::Vector{Int}
    num_layers::Int
    function GraphHist(a::Assignment{T,M}) where {T,M}
        θ = a.estimated_theta
        node_labels = a.node_labels
        new{typeof(θ[1]),M}(θ, node_labels, a.number_layers)
    end
end

"""
Network Histogram approximation [1].

Contains the estimated network histogram and the node labels.

# Fields
- `θ::Matrix{T}`: Estimated stochastic block model parameters.
- `node_labels::Vector{Int}`: Node labels for each node in the adjacency matrix used
    to estimate the network histogram.

# References
[1] - Olhede, Sofia C., and Patrick J. Wolfe. "Network histograms and universality of
blockmodel approximation." Proceedings of the National Academy of Sciences 111.41 (2014): 14722-14727.
"""
GraphHist


function get_moment_representation(g::GraphHist{T,2}) where {T}
    return g.θ
end



function get_moment_representation(g::GraphHist{T,3}) where {T}
    moments = zeros(size(g.θ,1), size(g.θ,2), 2^g.num_layers-1)
    transition = collect(kronecker([1 1; 0 1],g.num_layers))
    println(transition)
    for i in 1:size(g.θ,1)
        for j in 1:size(g.θ,2)
            moments[i,j,:] .= (transition*g.θ[i,j,:])[2:end]
        end
    end
    indices_for_moments = [findall(x -> x == 1, _index_to_binary(e, g.num_layers)) for e in 2:size(g.θ,3)]
    return moments, indices_for_moments
end
