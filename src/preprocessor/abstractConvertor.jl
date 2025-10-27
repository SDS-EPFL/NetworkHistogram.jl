abstract type AbstractConvertor end

Base.broadcastable(o::T) where {T <: AbstractConvertor} = Ref(o)

"""
    Convert data from its original form to a processed form suitable for SBM estimation.
"""
function (c::AbstractConvertor)(A; kwargs...)
    @error "to be implemented"
end

function to_distribution(c::AbstractConvertor, ps; kwargs...)
    @error "to be implemented"
end

include("categorical.jl")
include("continuous.jl")
