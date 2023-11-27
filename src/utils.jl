function laplacian(A)
    s = sum(A; dims = 1)
    return diagm(vec(s)) - A
end

function normalized_laplacian(A)
    L = zeros(size(A))
    degrees = vec(sum(A, dims = 1))
    for j in 1:size(A, 1)
        for i in 1:size(A, 2)
            if i == j
                L[i, j] = 1
            else
                L[i, j] = A[i, j] / sqrt(degrees[i] * degrees[j])
            end
        end
    end
    return L
end

function normalized_laplacian(A::AbstractArray{T, 3}) where {T}
    L = zeros(size(A, 1), size(A, 2))
    for layer in eachslice(A, dims = 3)
        L .+= normalized_laplacian(layer)
    end
    return L ./ size(A, 3)
end

function drop_disconnected_components(A::AbstractArray{T, 2}) where {T}
    indices = findall(x -> x != 0, vec(sum(A, dims = 1)))
    return A[indices, indices]
end

function drop_disconnected_components(A::AbstractArray{T, 3}) where {T}
    indices = findall(x -> x != 0, vec(sum(A, dims = (1, 3))))
    return A[indices, indices, :]
end

"""
    hamming_distance(x, y)

Compute the normalized Hamming distance between two vectors `x` and `y`.
"""
function hamming_distance(x::Vector{T}, y::Vector{T}) where {T}
    return sum(x .!= y) / length(x)
end

"""
    pairwise_hamming_distance(A)

Compute the pairwise Hamming distance between all rows of `A`. If `A` is a 3D
array, then the average Hamming distance for each layer of the array is returned.
"""
pairwise_hamming_distance

function pairwise_hamming_distance(matrix::AbstractArray{T, 2}) where {T}
    n = size(matrix, 1)
    dist_matrix = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        dist_matrix[i, j] = hamming_distance(matrix[i, :], matrix[j, :])
        dist_matrix[j, i] = dist_matrix[i, j]  # Symmetric matrix
    end
    return dist_matrix
end

function pairwise_hamming_distance(matrix::AbstractArray{T, 3}) where {T}
    n = size(matrix, 1)
    dist_matrix = zeros(n, n)
    for layer in eachslice(matrix, dims = 3)
        dist_matrix .+= pairwise_hamming_distance(layer)
    end
    return dist_matrix ./ size(matrix, 3)
end

function spectral_clustering(A, h)
    n = size(A, 1)

    L = 1 .- pairwise_hamming_distance(A) ./ n

    # Compute the degree matrix
    d = sum(L, dims = 2)

    # Compute the normalized Laplacian
    normalized_L = sum(1.0 ./ d) .* L .- sum(d) / sqrt(sum(d .^ 2))

    # Compute eigenvalues and eigenvectors of the normalized Laplacian

    decomp, history = partialschur(normalized_L, nev = 2, which = LR())
    _, eigen_vecs = partialeigen(decomp)

    # Extract the second eigenvector
    u = real.(eigen_vecs[:, 1])
    u = u .* sign(u[1]) # Set the first coordinate >= 0 wlog

    # Sort based on the embedding
    ind = sortperm(u, alg = QuickSort, rev = false)

    # Determine the number of clusters
    k = ceil(Int, n / h)

    # Initialize cluster assignments
    idxInit = zeros(Int, n)
    for i in 1:k
        for j in ((i - 1) * h + 1):min(n, i * h)
            idxInit[ind[j]] = i
        end
    end
    return idxInit
end
