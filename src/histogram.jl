struct GraphHist{T}
    θ::Matrix{T}
    node_labels::Vector{Int}
    function GraphHist(a::Assignment)
        θ = a.estimated_theta
        node_labels = a.node_labels
        new{typeof(θ[1])}(θ, node_labels)
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
