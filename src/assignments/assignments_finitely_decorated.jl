mutable struct AdditionalDataCategorical <: AdditionalData
    counts :: Matrix{Int}
    realized :: Array{Float64, 3}
    estimated_theta :: Array{Float64, 3}
    loglikelihood :: Float64
end


struct DefaultAssignmentCategorical <: Assignment
    group_size :: GroupSize{Int}
    node_labels :: Vector{Int}
    additional_data :: AdditionalDataCategorical
end
