abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end

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
objecta vector of node labels.

# Implemenented rules
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

# check https://github.com/TrainOfCode/LocalFennelPartitioning.jl/tree/main
# check https://github.com/JuliaSparse/Metis.jl
