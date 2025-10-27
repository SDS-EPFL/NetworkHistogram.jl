### =======================================================================================
###                         Categorical Convertor
### =======================================================================================

struct CategoricalConvertor{T} <: AbstractConvertor
    m::Int  # number of categories
    map::Dict{T, Int}
end

function CategoricalConvertor(data::AbstractArray{T}) where {T}
    categories = sort(unique(data))
    m = length(categories)
    map = Dict{T, Int}(categories[i] => i for i in 1:m)
    return CategoricalConvertor{T}(m, map)
end

function num_bins(c::CategoricalConvertor)
    return c.m
end

function (c::CategoricalConvertor)(obs::T) where {T}
    return c.map[obs]
end

function to_distribution(
        c::CategoricalConvertor{T}, ps::AbstractVector{T2}; kwargs...) where {T, T2}
    @argcheck length(ps)==c.m "Length of probabilities must match number of categories"
    support = sort(collect(keys(c.map)))
    probabilities = SVector{c.m, T2}(ps[c.map[s]] for s in support)
    return DiscreteNonParametric(support, probabilities)
end
