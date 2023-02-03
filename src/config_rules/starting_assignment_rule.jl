abstract type StartingAssignment end
struct OrderedStart <: StartingAssignment end
struct RandomStart <: StartingAssignment end

function initialise_node_labels(A, h, ::OrderedStart)
    group_size = GroupSize(size(A,1), h)
    inverse_rle(1:length(group_size), group_size)
    return node_labels, group_size
end

function initialise_node_labels(A, h, ::RandomStart)
    node_labels, group_size = initialise_node_labels(A, h, OrderedStart())
    permute!(node_labels, rand(1:length(node_labels)))
    return node_labels, group_size
end