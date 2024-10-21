abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end
struct SpectralStart <: StartingAssignment end
struct MetisStart <: StartingAssignment end
struct FromAssignment{A} <: StartingAssignment
    assignment::A
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
        Int, Metis.partition(Metis.graph(g.graph), length(group_size)))
    check_compatiblity(group_size, node_labels)
    return group_size, node_labels
end

function initialize_node_labels(g, h, rule::FromAssignment{A}) where {A}
    group_size = GroupSize(number_nodes(g), h)
    check_compatiblity(group_size, rule.assignment.node_labels)
    return group_size, rule.assignment.node_labels
end
