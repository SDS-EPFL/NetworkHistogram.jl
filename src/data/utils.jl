function drop_isolated_vertices(A)
    degrees = vec(sum(A, dims = 2))
    return A[degrees .> 0, degrees .> 0]
end

####### From Erdos.jl #######

function minutype(n::Integer)
    @assert n ≥ 0
    if n < 2^8
        return UInt8
    elseif n < 2^16
        return UInt16
    elseif n < 2^32
        return UInt32
    elseif n < 2^64
        return UInt64
    end
    error("No type big enough")
end
