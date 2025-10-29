struct BinaryConvertor <: AbstractConvertor end

function (c::BinaryConvertor)(obs::T) where {T <: Union{Real, Bool}}
    return obs == 1 ? true : false
end

function to_distribution(
        c::BinaryConvertor, p::T; kwargs...) where {T <: Real}
    return p
end
