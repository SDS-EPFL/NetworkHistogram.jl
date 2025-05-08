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
    apply_swap!(assignment, swap)
    assignment.θ = deepcopy(swap.workspace.θ)
    assignment.log_likelihood = deepcopy(swap.workspace.log_likelihood_per_group)
end

function swap_node_labels!(a::Assignment, i, j)
    a.node_labels[i], a.node_labels[j] = a.node_labels[j], a.node_labels[i]
end

function apply_swap!(a::Assignment, s::Swap)
    g1 = group(a, s.u)
    g2 = group(a, s.v)
    groups_concerned = Set([minmax(g1, g2)])
    for (u, g_old, g_new) in [(s.u, g1, g2), (s.v, g2, g1)]
        # iterate over neighbors of u and get the decoration of the edge
        for (v,d) in iterate_neighbors(a.dists, u)
            g_v = group(a, v)
            a.θ[g_old, g_v] = remove_from(a.θ[g_old, g_v], d)
            a.θ[g_new, g_v] = add_to(a.θ[g_new, g_v], d)
            push!(groups_concerned, minmax(g_new, g_v))
            push!(groups_concerned, minmax(g_old, g_v))
        end
    end
    @show a.θ
    fast_ll_update!(a, groups_concerned)

    swap_node_labels!(a, s.u, s.v)
end


## below can be specialised for Bernoulli probably

function fast_ll_update!(a, groups_concerned)
    for g in groups_concerned
        # Use get_edges_in_groups to get the correct set of edges
        edges = get_edges_in_groups(a, g[1], g[2])
        a.log_likelihood[g[1], g[2]] = loglikelihood(a.θ[g[1], g[2]], edges)
    end
end
