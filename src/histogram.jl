struct GraphHist{T} <: AbstractMatrix{T}
    θ::Matrix{T}
    node_labels::Vector{Int}
    function GraphHist(a::Assignment)
        θ = a.estimated_theta
        node_labels = a.node_labels
        new{typeof(θ[1])}(θ, node_labels)
    end
end
