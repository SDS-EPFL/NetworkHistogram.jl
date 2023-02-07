function laplacian(A)
    D = diagm(vec(sum(A, dims = 2)))
    return D - A
end

function normalized_laplacian(A)
    @warn "Suppose no isolated vertices, no check were performed before"
    D_inv_half = diagm(vec(1 ./ sqrt.(sum(A, dims = 2))))
    return I - D_inv_half * A * D_inv_half
end

function spectral_clustering(A, K)
    L = normalized_laplacian(A)
    eigenvalues, eigenvectors = eigen(L)
    eigenvectors = eigenvectors[:, 2:(K + 1)]
    node_labels = kmeans(eigenvectors, K).assignments
    group_size = GroupSize(size(A, 1), K)
    return node_labels, group_size
end
