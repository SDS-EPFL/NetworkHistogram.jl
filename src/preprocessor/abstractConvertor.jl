abstract type AbstractConvertor end

Base.broadcastable(o::AbstractConvertor) = Ref(o)

"""
    Convert data from its original form to a processed form suitable for SBM estimation.
"""
function (c::AbstractConvertor)(A; kwargs...)
    @error "to be implemented"
end

function to_distribution(c::AbstractConvertor, ps; kwargs...)
    @error "to be implemented"
end

get_convertor(s::String, ; kwargs...) = get_convertor(Symbol(s); kwargs...)
get_convertor(s::Symbol; kwargs...) = get_convertor(Val(s); kwargs...)
get_convertor(::T; kwargs...) where {T} = @error "No convertor found for type $T"

include("categorical.jl")
include("continuous.jl")

function get_convertor(::Val{:categorical}; kwargs...)
    return CategoricalConvertor(kwargs[:num_categories])
end

function get_convertor(::Val{:continuous}; kwargs...)
    return UnitIntervalConvertor(kwargs[:num_bins])
end
