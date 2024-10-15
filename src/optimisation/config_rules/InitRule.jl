abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end


struct InitRule{S <: StartingAssignment, I}
    starting_assignment_rule::S
    assignment_rule::I
end

function make_assignment(A, h, init_rule::InitRule{S, Nothing}) where S
    return Assignment(initialize_node_labels(A, h, init_rule.starting_assignment_rule)...)
end

"""
    initialize_node_labels(A, h, starting_assignment_rule::StartingAssignment)

initialize node labels based on the `starting_assignment_rule`, and return a vector of
node labels and a `GroupSize` object.

# Implemenented rules
- `OrderedStart()`: Sequentially assign nodes to groups based on the ordering of `A`.
- `RandomStart()`: Randomly assign nodes to groups.
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

# check https://github.com/TrainOfCode/LocalFennelPartitioning.jl/tree/main
# check https://github.com/JuliaSparse/Metis.jl
