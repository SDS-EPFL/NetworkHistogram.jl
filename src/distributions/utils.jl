const logtwo = log(2.0)

sumlog(x::AbstractArray{<:Real}) = sum(log,x)

function sumlog(x::AbstractArray{<:AbstractFloat})
    sig = one(T)
    ex = zero(exponent(one(T)))
    bound = floatmax(T) / 2
    for xj in x
        sig *= significand(xj)
        ex += exponent(xj)
        if sig > bound
            (a, b) = (significand(sig), exponent(sig))
            sig = a
            ex += b
        end
    end
    log(sig) + logtwo * ex
end
