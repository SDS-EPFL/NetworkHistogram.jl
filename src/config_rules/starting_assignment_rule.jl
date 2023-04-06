abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end
struct EigenStart <: StartingAssignment end


"""
    initialise_node_labels(A, h, starting_assignment_rule::StartingAssignment)

Initialise node labels based on the `starting_assignment_rule`, and return a vector of
node labels and a `GroupSize` object.

# Implemenented rules
- `OrderedStart()`: Sequentially assign nodes to groups based on the ordering of `A`.
- `RandomStart()`: Randomly assign nodes to groups.
- `EigenStart()`: Assign nodes to groups based on the second eigenvector of the normalized Laplacian.
"""
initialise_node_labels


function initialise_node_labels(A, h, ::OrderedStart)
    group_size = GroupSize(size(A, 1), h)
    node_labels = inverse_rle(1:length(group_size), group_size)
    return node_labels, group_size
end

function initialise_node_labels(A, h, ::RandomStart)
    node_labels, group_size = initialise_node_labels(A, h, OrderedStart())
    node_labels = shuffle!(node_labels)
    return node_labels, group_size
end


function initialise_node_labels(A, h, ::EigenStart)
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
