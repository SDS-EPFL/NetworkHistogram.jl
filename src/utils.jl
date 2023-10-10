function laplacian(A)
    s = sum(A; dims = 1)
    return diagm(vec(s)) - A
end

function normalized_laplacian(A)
    D_scale = vec(sum(A; dims = 1)) .^ -0.5
    return I(size(A, 1)) - Diagonal(D_scale) * (A * Diagonal(D_scale))
end

function normalized_laplacian(A::Array{T, 3}) where {T}
    D_scale = vec(sum(A; dims = [1, 3])) .^ -0.5
    return I(size(A, 1)) -
           Diagonal(D_scale) * (dropdims(sum(A, dims = 3), dims = 3) * Diagonal(D_scale))
end
