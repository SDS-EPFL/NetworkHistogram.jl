mutable struct WorkspaceSwap{D,F}
    θ::SymArray{D}
    log_likelihood_per_group::SymArray{F}
end

mutable struct Swap{W}
    u::Int
    v::Int
    workspace::W
end


function make_swap(a::Assignment, id)
    return Swap(id[1], id[2], WorkspaceSwap(deepcopy(a.θ), deepcopy(a.log_likelihood)))
end

function make_swap!(swap::Swap, a::Assignment, id)
    swap.u, swap.v = id
    swap.workspace.θ = deepcopy(a.θ)
    swap.workspace.log_likelihood_per_group = deepcopy(a.log_likelihood)
end

function revert_swap!(assignment::Assignment, swap::Swap)
    # swap labels back to original
    swap_node_labels!(assignment, swap.u, swap.v)
    # restore saved θ and log likelihoods
    assignment.θ = deepcopy(swap.workspace.θ)
    assignment.log_likelihood = deepcopy(swap.workspace.log_likelihood_per_group)
end

function swap_node_labels!(a::Assignment, i, j)
    a.node_labels[i], a.node_labels[j] = a.node_labels[j], a.node_labels[i]
end

# for reference and testing
function _slow_swap!(a::Assignment, s::Swap)
    swap_node_labels!(a, s.u, s.v)
    a.θ, a.log_likelihood = _compute_theta_and_ll(a.node_labels, a.dists, a.edges, a.θ[1,1])
end

# apply_swap!(a::Assignment, s::Swap) = _slow_swap!(a, s)


function apply_swap!(a::Assignment, s::Swap)
    u,v = s.u, s.v
    gu = a.node_labels[u]
    gv = a.node_labels[v]
    groups_concerned = Set{Tuple{Int,Int}}([minmax(gu, gv)])

    for (node, d) in iterate_neighbors(a.dists, u)
        if node == v
            continue
        end
        g1 = a.node_labels[node]
        a.θ[gv, g1] = add_to(a.θ[gv, g1], d)
        a.θ[gu, g1] = remove_from(a.θ[gu, g1], d)
        push!(groups_concerned, minmax(gu,g1))
        push!(groups_concerned, minmax(gv,g1))
    end

    for (index, (node, d)) in enumerate(iterate_neighbors(a.dists, v))
        if node == u
            continue
        end
        g2 = a.node_labels[node]
        a.θ[gu, g2] = add_to(a.θ[gu, g2], d)
        a.θ[gv, g2] = remove_from(a.θ[gv, g2], d)
        push!(groups_concerned, minmax(gv,g2))
        push!(groups_concerned, minmax(gu,g2))
    end

    swap_node_labels!(a, u, v)
    for (g1, g2) in groups_concerned
        a.log_likelihood[g1, g2] = 0.0
        for e in get_edges_in_groups(a.node_labels, a.edges, g1, g2)
            a.log_likelihood[g1, g2] += logpdf(a.θ[g1, g2], e)
        end
    end
end
