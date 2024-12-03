using CategoricalArrays
using CategoricalDistributions
using Distributions


struct Encoder{F, S}
    breaks::Vector{F}
    labels::Vector{S}
    extended::Bool

    function Encoder(
            _breaks::AbstractVector{F}, labels = CategoricalArrays.default_formatter;
            extend = missing) where {F}
        breaks = sort(_breaks)
        n = length(breaks)
        from = breaks[1:(n - 1)]
        to = breaks[2:n]
        firstlevel = labels(from[1], to[1], 1,
            leftclosed = breaks[1] != breaks[2], rightclosed = false)
        levs = Vector{typeof(firstlevel)}(undef, n - 1)
        levs[1] = firstlevel
        for i in 2:(n - 2)
            levs[i] = labels(from[i], to[i], i,
                leftclosed = breaks[i] != breaks[i + 1], rightclosed = false)
        end
        levs[end] = labels(from[end], to[end], n - 1,
            leftclosed = breaks[end - 1] != breaks[end],
            rightclosed = coalesce(extend, false))

        new{F, typeof(firstlevel)}(breaks, levs, coalesce(extend, false))
    end
end

function convert(encoder::Encoder, x::T) where {T <: Real}
    if x < encoder.breaks[1] || x > encoder.breaks[end]
        throw(ArgumentError("Value $x out of bounds $(encoder.breaks[1]) - $(encoder.breaks[end])"))
    end
    if x == encoder.breaks[end] && encoder.extended
        return encoder.labels[end]
    end
    if x == encoder.breaks[1]
        return encoder.labels[1]
    end
    return encoder.labels[findlast(y-> y <= x, encoder.breaks)]
end


function convert(encoder::Encoder, x::String)
    index = findfirst(l -> l == x, encoder.labels)
    if isnothing(index)
        throw(ArgumentError("Value $x not found in $(encoder.labels)"))
    end
    if index == 1
        return encoder.breaks[1], encoder.breaks[2]
    elseif index == length(encoder.labels)
        return encoder.breaks[end-1], encoder.breaks[end]
    else
        return encoder.breaks[index-1], encoder.breaks[index]
    end
end


struct DiscretisedDist{S, F, L}
    dist::UnivariateFinite{S}
    encoding::Encoder{F,L}
end
