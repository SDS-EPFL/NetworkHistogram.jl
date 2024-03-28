abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end
struct EigenStart <: StartingAssignment end
struct DistStart <: StartingAssignment end
struct LSBM <: StartingAssignment end

"""
    initialize_node_labels(A, h, starting_assignment_rule::StartingAssignment)

initialize node labels based on the `starting_assignment_rule`, and return a vector of
node labels and a `GroupSize` object.

# Implemenented rules
- `OrderedStart()`: Sequentially assign nodes to groups based on the ordering of `A`.
- `RandomStart()`: Randomly assign nodes to groups.
- `EigenStart()`: Assign nodes to groups based on the second eigenvector of the normalized Laplacian.
- `DistStart()`: Assign nodes to groups based on the Hamming distance between rows of `A`.
"""
initialize_node_labels

function initialize_node_labels(A, h, ::OrderedStart)
    group_size = GroupSize(size(A, 1), h)
    node_labels = inverse_rle(1:length(group_size), group_size)
    return node_labels, group_size
end

function initialize_node_labels(A, h, ::RandomStart)
    node_labels, group_size = initialize_node_labels(A, h, OrderedStart())
    node_labels = shuffle!(node_labels)
    return node_labels, group_size
end

function initialize_node_labels(A, h, ::EigenStart)
    group_size = GroupSize(size(A, 1), h)
    node_labels = zeros(Int, size(A, 1))

    laplacian = normalized_laplacian(A)
    _, eigenvectors = eigs(laplacian, nev = 2, which = :LR, tol = 1e-2)
    #_, eigenvectors = eigen(Symmetric(laplacian), (size(A, 1) - 1):(size(A, 1) - 1))

    # get 2nd eigenvector, sort its components
    indices = sortperm(eigenvectors[:, 1])
    # bin them into groups of correct size
    start = 1
    for (i, group) in enumerate(group_size)
        stop = start + group - 1
        node_labels[indices[start:stop]] .= i
        start = stop + 1
    end
    return node_labels, group_size
end

function initialize_node_labels(A, h, ::DistStart)
    group_size = GroupSize(size(A, 1), h)
    node_labels = spectral_clustering(A, h)
    return node_labels, group_size
end


function initialize_node_labels(A::Array{I, 3}, h, ::LSBM; r=2) where {I}
    @error "LSBM starting point does not guarantee equally sized groups, balancing is to be implemented"
    # init containers
    group_size = GroupSize(size(A, 1), h)
    node_labels = zeros(Int, size(A, 1))

    # random weight function with 0 at the first index (all 0 vector)
    random_function = rand(2^size(A, 3))
    random_function[1] = 0

    # build weighted adjacency matrix
    A_categorical = update_adj(A)
    weighted_adjacency = zeros(Float64, size(A, 1), size(A, 2))
    for j in 1:size(A,2)
        for i in 1:size(A,1)
            if i != j
                weighted_adjacency[j,i] = random_function[A_categorical[j,i]]
            end
        end
    end

    eigenvals, eigenvectors = eigs(weighted_adjacency, nev = r, which = :LR, tol = 1e-2)
    z = eigenvectors .* (eigenvals'/eigenvals[1])

    result = kmeans(transpose(z), size(group_size,1))
    node_labels = result.assignments
    return node_labels, group_size
end
