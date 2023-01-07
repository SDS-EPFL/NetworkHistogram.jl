struct Assignment
    number_nodes::Int
    number_groups::Int
    node_labels::Vector{Int}

    counts::Matrix{Int}
    realized::Matrix{Float64}
    estimated_theta::Matrix{Float64}

    likelihood::Float64
    function Assignment(A, h)

        # need to fill this in

        new(number_nodes,
            number_groups,
            node_labels,
            counts,
            realized,
            estimated_theta,
            likelihood)
    end

    function Assignment(a::Assignment, likelihood::Real)
        new(a.number_nodes,
            a.number_groups,
            a.node_labels,
            a.counts,
            a.realized,
            a.estimated_theta,
            likelihood)
    end
end

function initialize(A::Matrix{Int}, h::Int)
    old_store = Assignment(A, h)
    new_store = deepcopy(oldstore)
    history = MVHistory([
                            :likelihood => QHistory(Float64),
                            :best_likelihood => QHistory(Float64),
                            :proposal_likelihood => QHistory(Float64),
                        ])
    return old_store, new_store, history
end
