struct Assignment
    number_nodes::Int
    number_groups::Int
    node_labels::Vector{Int}

    counts::Matrix{Int}
    realized::Matrix{Float64}
    estimated_theta::Matrix{Float64}

    likelihood::Float64
    function Assignment(A, node_labels)
        
        number_groups = length(unique(node_labels))
        number_nodes = length(node_labels)
        estimated_theta = zeros(Float64, number_nodes, number_groups)
        counts = zeros(Int64, number_groups, number_groups)
        realized = zeros(Int64, number_groups, number_groups)
        size_groups = [count(==(element), node_labels) for element in unique(node_labels)]

        @inbounds @simd for k in 1:number_groups
            for l in k:number_groups
                realized[k, l] = sum(adjacency_matrix[node_labels.==k, node_labels.==l])
                realized[l, k] = realized[k, l]
                counts[k, l] = size_groups[k] * size_groups[l]
                counts[l, k] = counts[k, l]
            end
        end

        @inbounds @simd for k in 1:number_groups
            counts[k, k] = size_groups[k] * (size_groups[k] - 1) / 2
            realized[k, k] = sum(adjacency_matrix[node_labels.==k, node_labels.==k]) / 2
        end

        estimated_theta = realized ./ counts
      
        new(
            number_nodes,
            number_groups,
            node_labels,
            counts,
            realized,
            estimated_theta,
            likelihood
        )
    end

    function Assignment(a::Assignment, likelihood)
        new(
            a.number_nodes,
            a.number_groups,
            a.node_labels,
            a.counts,
            a.realized,
            a.estimated_theta,
            likelihood
        )
    end

end

function initialize(A, h)
    node_labels = initialise_node_labels(A, h)
    old_store = Assignment(A, node_labels)
    new_store = deepcopy(oldstore)
    history = MVHistory([:likelihood => QHistory(Float64), :best_likelihood => QHistory(Float64), :proposal_likelihood => QHistory(Float64)])
    return old_store, new_store, history
end
