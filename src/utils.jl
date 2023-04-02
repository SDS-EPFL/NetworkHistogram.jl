function laplacian(A)
    s = sum(A; dims = 1)
    return diagm(vec(s)) - A
end

function normalized_laplacian(A)
    D_scale = vec(sum(A; dims = 1)) .^ -0.5
    return I(size(A, 1)) - Diagonal(D_scale) * (A * Diagonal(D_scale))
end
