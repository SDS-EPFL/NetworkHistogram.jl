mutable struct AdditionalDataBinary <: AdditionalDataBinary
    counts :: Matrix{Int}
    realized :: Matrix{Float64}
    estimated_theta :: Matrix{Float64}
    loglikelihood::Float64
end

struct DefaultAssignmentBinary <: AssignmentBinary
    group_size :: GroupSize{Int}
    node_labels :: Vector{Int}
    additional_data :: AdditionalDataBinary
end
