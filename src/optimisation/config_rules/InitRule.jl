abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end
struct SpectralStart <: StartingAssignment end
struct MetisStart <: StartingAssignment end
struct FromAssignment{A} <: StartingAssignment
    assignment::A
end
struct HigherOrderSpectralStart <: StartingAssignment
    k::Int
end

struct InitRule{S <: StartingAssignment, I}
    starting_assignment_rule::S
    assignment_rule::I
end

function make_assignment(g, h, init_rule::InitRule{S, Nothing}) where {S}
    return Assignment(initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)...)
end

"""
    initialize_node_labels(g, h, starting_assignment_rule::StartingAssignment)

initialize node labels based on the `starting_assignment_rule`, and return a `GroupSize`
object and a vector of node labels.

# Implemented rules
- `OrderedStart()`: Sequentially assign nodes to groups based on the ordering of `A`.
- `RandomStart()`: Randomly assign nodes to groups.
- `SpectralStart()`: Assign nodes to groups based on spectral clustering.
- `MetisStart()`: Assign nodes to groups based on Metis partitioning.
- `FromAssignment(a)`: Assign nodes to groups based on the given assignment `a`.
"""
initialize_node_labels

function initialize_node_labels(g, h, ::OrderedStart)
    group_size = GroupSize(number_nodes(g), h)
    node_labels = StatsBase.inverse_rle(1:length(group_size), group_size)
    return group_size, node_labels
end

function initialize_node_labels(g, h, ::RandomStart)
    group_size, node_labels = initialize_node_labels(g, h, OrderedStart())
    Random.shuffle!(node_labels)
    return group_size, node_labels
end

function initialize_node_labels(g, h, ::SpectralStart)
    group_size = GroupSize(number_nodes(g), h)
    node_labels = zeros(Int, number_nodes(g))

    laplacian = normalized_laplacian(g)
    _, eigenvectors = Arpack.eigs(laplacian, nev = 2, which = :LR)
    # get 2nd eigenvector, sort its components
    indices = sortperm(eigenvectors[:, 1])
    # bin them into groups of correct size
    start = 1
    for (i, group) in enumerate(group_size)
        stop = start + group - 1
        node_labels[indices[start:stop]] .= i
        start = stop + 1
    end
    return group_size, node_labels
end

function initialize_node_labels(g, h, ::MetisStart)
    group_size = GroupSize(number_nodes(g), h)
    node_labels = convert.(
        Int, Metis.partition(Metis.graph(g), length(group_size)))
    check_compatiblity!(node_labels, group_size)
    return group_size, node_labels
end

function initialize_node_labels(g, h, rule::FromAssignment{A}) where {A}
    group_size = GroupSize(number_nodes(g), h)
    check_compatiblity!(rule.assignment.node_labels, group_size)
    return group_size, rule.assignment.node_labels
end

function initialize_node_labels(g, h, rule::HigherOrderSpectralStart)
    throw(ArgumentError("Not implemented yet, need to finish with Clustering.jl"))
    group_size = GroupSize(number_nodes(g), h)
    laplacian = normalized_laplacian(g)
    results = IterativeSolvers.lobpcg(laplacian, true, rule.k)
    return group_size, node_labels
end
