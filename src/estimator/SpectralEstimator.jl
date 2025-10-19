"""
    SpectralEstimator{T}

Spectral clustering estimator for Stochastic Block Models (SBM).

This estimator uses spectral clustering to partition nodes into groups based on
the graph structure. It computes the normalized Laplacian and performs k-means
clustering on the top k eigenvectors.

# Fields
- `k::Int`: Number of groups/communities to detect
- `eig_size::Int`: Number of eigenvectors to use (default: k)
- `adjacency_type::Symbol`: Type of adjacency matrix to use (`:binary`, `:weighted`)
- `laplacian_type::Symbol`: Type of Laplacian (`:normalized`, `:unnormalized`)
- `max_kmeans_iter::Int`: Maximum iterations for k-means clustering
- `balanced::Bool`: If true, forces balanced community sizes (default: false)

# Example
```julia
# Binary adjacency matrix with balanced communities
A = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]
estimator = SpectralEstimator(2, balanced=true)
labels = estimate(estimator, A)
```
"""
struct SpectralEstimator{T <: Real} <: SBMEstimator
    k::Int
    eig_size::Int
    adjacency_type::Symbol
    laplacian_type::Symbol
    max_kmeans_iter::Int
    balanced::Bool

    function SpectralEstimator(k::Int;
            eig_size::Int = k,
            adjacency_type::Symbol = :binary,
            laplacian_type::Symbol = :normalized,
            max_kmeans_iter::Int = 100,
            balanced::Bool = false)
        @argcheck k>0 "Number of groups k must be positive"
        @argcheck adjacency_type in [:binary, :weighted] "adjacency_type must be :binary or :weighted"
        @argcheck laplacian_type in [:normalized, :unnormalized] "laplacian_type must be :normalized or :unnormalized"
        @argcheck max_kmeans_iter>0 "max_kmeans_iter must be positive"
        new{Float64}(k, eig_size, adjacency_type, laplacian_type, max_kmeans_iter, balanced)
    end
end

"""
    estimate(estimator::SpectralEstimator, data; progress = true)

Perform spectral clustering on the network data.

# Arguments
- `estimator::SpectralEstimator`: The spectral estimator configuration
- `data`: The adjacency matrix or network data
- `progress::Bool`: Whether to show progress information (for compatibility with other estimators)

# Returns
- `labels::Vector{Int}`: Node group assignments (1 to k)

# Algorithm
1. Construct adjacency matrix from data
2. Compute the Laplacian matrix (normalized or unnormalized)
3. Compute eigenvectors corresponding to smallest eigenvalues
4. Perform k-means clustering on the eigenvectors
5. Return cluster assignments
"""
function estimate(estimator::SpectralEstimator, data; progress = true)
    progress && @info "Starting spectral clustering with k=$(estimator.k)"

    # Convert data to adjacency matrix
    A = construct_adjacency(data, estimator.adjacency_type)
    n = size(A, 1)

    @argcheck n>=estimator.k "Number of nodes ($n) must be >= number of groups ($(estimator.k))"

    # Compute Laplacian
    L = compute_laplacian(A, estimator.laplacian_type)

    # Compute eigenvectors
    progress && @info "Computing eigenvectors..."
    eigvals, eigvecs = compute_spectral_embedding(L, estimator.eig_size)

    # Normalize rows for normalized spectral clustering
    if estimator.laplacian_type == :normalized
        eigvecs = normalize_rows(eigvecs)
    end

    # Perform k-means clustering
    progress && @info "Performing k-means clustering..."
    if estimator.balanced
        labels = balanced_kmeans_clustering(eigvecs, estimator.k, estimator.max_kmeans_iter)
    else
        labels = kmeans_clustering(eigvecs, estimator.k, estimator.max_kmeans_iter)
    end

    progress && @info "Spectral clustering complete"
    return labels
end

"""
    estimate(estimator::SpectralEstimator, data, initial_labels; progress = true)

Perform spectral clustering on the network data. The initial_labels are ignored
as spectral clustering doesn't use an initialization.

This method signature is provided for compatibility with other estimators.
"""
function estimate(estimator::SpectralEstimator, data, initial_labels; progress = true)
    return estimate(estimator, data; progress = progress)
end

"""
    construct_adjacency(data, adjacency_type::Symbol)

Construct an adjacency matrix from the data.

# Arguments
- `data`: Network data (can be a matrix with various edge types)
- `adjacency_type`: Either `:binary` or `:weighted`

# Returns
- Symmetric adjacency matrix
"""
function construct_adjacency(data::AbstractMatrix, adjacency_type::Symbol)
    n = size(data, 1)
    A = zeros(Float64, n, n)

    if adjacency_type == :binary
        # Binary adjacency: edge exists if data is not nothing/zero
        for i in 1:n
            for j in (i + 1):n
                if !isnothing(data[i, j]) && data[i, j] != 0
                    A[i, j] = 1.0
                    A[j, i] = 1.0
                end
            end
        end
    elseif adjacency_type == :weighted
        # Weighted adjacency: use the actual values
        for i in 1:n
            for j in (i + 1):n
                if !isnothing(data[i, j])
                    if data[i, j] isa AbstractArray
                        # For categorical data, use sum or count
                        weight = sum(data[i, j])
                    else
                        weight = float(data[i, j])
                    end
                    A[i, j] = weight
                    A[j, i] = weight
                end
            end
        end
    end

    return A
end

"""
    compute_laplacian(A::AbstractMatrix, laplacian_type::Symbol)

Compute the graph Laplacian matrix.

# Arguments
- `A`: Adjacency matrix
- `laplacian_type`: Either `:normalized` or `:unnormalized`

# Returns
- Laplacian matrix
"""
function compute_laplacian(A::AbstractMatrix, laplacian_type::Symbol)
    n = size(A, 1)
    d = vec(sum(A, dims = 2))  # Degree vector

    if laplacian_type == :unnormalized
        # L = D - A
        D = Diagonal(d)
        return D - A
    elseif laplacian_type == :normalized
        # L = I - D^{-1/2} A D^{-1/2}
        # Handle zero degrees
        d_inv_sqrt = zeros(n)
        for i in 1:n
            d_inv_sqrt[i] = d[i] > 0 ? 1.0 / sqrt(d[i]) : 0.0
        end
        D_inv_sqrt = Diagonal(d_inv_sqrt)
        return I - D_inv_sqrt * A * D_inv_sqrt
    end
end

"""
    compute_spectral_embedding(L::AbstractMatrix, k::Int)

Compute the spectral embedding by finding eigenvectors corresponding to
the k smallest eigenvalues of the Laplacian.

# Arguments
- `L`: Laplacian matrix
- `k`: Number of eigenvectors to compute

# Returns
- `eigvals`: The k smallest eigenvalues
- `eigvecs`: Matrix where each row is a node and columns are eigenvector components
"""
function compute_spectral_embedding(L::AbstractMatrix, k::Int)
    # Compute smallest k eigenvalues and eigenvectors
    # Use eigen for small matrices, could use iterative methods for large ones
    n = size(L, 1)

    if n <= 1000
        # For small matrices, compute all eigenvalues
        F = eigen(Symmetric(L))
        idx = sortperm(F.values)[1:k]
        return F.values[idx], F.vectors[:, idx]
    else
        # For larger matrices, use iterative solver (if available)
        # For now, still use full eigen but this could be optimized
        F = eigen(Symmetric(L))
        idx = sortperm(F.values)[1:k]
        return F.values[idx], F.vectors[:, idx]
    end
end

"""
    normalize_rows(X::AbstractMatrix)

Normalize each row of the matrix to unit length.

# Arguments
- `X`: Matrix to normalize

# Returns
- Matrix with normalized rows
"""
function normalize_rows(X::AbstractMatrix)
    n, k = size(X)
    X_norm = similar(X)

    for i in 1:n
        row_norm = norm(X[i, :])
        if row_norm > 0
            X_norm[i, :] = X[i, :] / row_norm
        else
            X_norm[i, :] = X[i, :]
        end
    end

    return X_norm
end

"""
    kmeans_clustering(X::AbstractMatrix, k::Int, max_iter::Int)

Perform k-means clustering on the rows of X.

# Arguments
- `X`: Data matrix (n × d), where n is number of points, d is dimensionality
- `k`: Number of clusters
- `max_iter`: Maximum number of iterations

# Returns
- `labels::Vector{Int}`: Cluster assignments (1 to k)
"""
function kmeans_clustering(X::AbstractMatrix, k::Int, max_iter::Int)
    n, d = size(X)

    # Initialize centers by randomly selecting k rows
    center_indices = randperm(n)[1:k]
    centers = X[center_indices, :]
    labels = zeros(Int, n)

    converged = false
    for iter in 1:max_iter
        # Assignment step
        old_labels = copy(labels)
        for i in 1:n
            min_dist = Inf
            best_cluster = 1
            for j in 1:k
                dist = sum(abs2, X[i, :] - centers[j, :])
                if dist < min_dist
                    min_dist = dist
                    best_cluster = j
                end
            end
            labels[i] = best_cluster
        end

        # Check convergence
        if labels == old_labels
            converged = true
            break
        end

        # Update step
        for j in 1:k
            cluster_points = findall(labels .== j)
            if !isempty(cluster_points)
                centers[j, :] = vec(mean(X[cluster_points, :], dims = 1))
            end
        end
    end

    return labels
end

"""
    balanced_kmeans_clustering(X::AbstractMatrix, k::Int, max_iter::Int)

Perform balanced k-means clustering on the rows of X, ensuring approximately equal-sized clusters.

This uses a greedy assignment approach where each cluster is filled to its target size
by assigning the closest points to each cluster's center, respecting size constraints.

# Arguments
- `X`: Data matrix (n × d), where n is number of points, d is dimensionality
- `k`: Number of clusters
- `max_iter`: Maximum number of iterations

# Returns
- `labels::Vector{Int}`: Cluster assignments (1 to k) with balanced sizes
"""
function balanced_kmeans_clustering(X::AbstractMatrix, k::Int, max_iter::Int)
    n, d = size(X)
    target_size = n ÷ k
    remainder = n % k

    # Initialize centers by randomly selecting k rows
    center_indices = randperm(n)[1:k]
    centers = X[center_indices, :]
    labels = zeros(Int, n)

    for iter in 1:max_iter
        old_labels = copy(labels)

        # Compute all distances
        distances = zeros(n, k)
        for i in 1:n
            for j in 1:k
                distances[i, j] = sum(abs2, X[i, :] - centers[j, :])
            end
        end

        # Balanced assignment using greedy approach
        labels = balanced_assignment(distances, k, target_size, remainder)

        # Check convergence
        if labels == old_labels
            break
        end

        # Update centers
        for j in 1:k
            cluster_points = findall(labels .== j)
            if !isempty(cluster_points)
                centers[j, :] = vec(mean(X[cluster_points, :], dims = 1))
            end
        end
    end

    return labels
end

"""
    balanced_assignment(distances::Matrix, k::Int, target_size::Int, remainder::Int)

Assign points to clusters in a balanced way using a greedy approach.

Each cluster gets exactly `target_size` or `target_size + 1` points (depending on remainder).

# Arguments
- `distances`: Matrix of distances from each point to each cluster center (n × k)
- `k`: Number of clusters
- `target_size`: Base number of points per cluster
- `remainder`: Number of clusters that get one extra point

# Returns
- `labels::Vector{Int}`: Balanced cluster assignments
"""
function balanced_assignment(distances::Matrix, k::Int, target_size::Int, remainder::Int)
    n = size(distances, 1)
    labels = zeros(Int, n)
    cluster_sizes = zeros(Int, k)
    max_sizes = fill(target_size, k)
    max_sizes[1:remainder] .+= 1

    # Create a list of (distance, point_idx, cluster_idx) tuples
    assignments = []
    for i in 1:n
        for j in 1:k
            push!(assignments, (distances[i, j], i, j))
        end
    end

    # Sort by distance (greedy: assign closest points first)
    sort!(assignments, by = x -> x[1])

    # Assign points greedily while respecting size constraints
    assigned = falses(n)
    for (dist, i, j) in assignments
        if !assigned[i] && cluster_sizes[j] < max_sizes[j]
            labels[i] = j
            cluster_sizes[j] += 1
            assigned[i] = true
        end
        if all(assigned)
            break
        end
    end

    return labels
end
