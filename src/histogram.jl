struct GraphHist <: AbstractMatrix{T}
    θ::Matrix{Float64}
    node_labels::Vector{Int}
    function GraphHist(a::Assignment)
        θ = a.estimated_theta
        node_labels = a.node_labels
        new(θ,node_labels)
    end
end