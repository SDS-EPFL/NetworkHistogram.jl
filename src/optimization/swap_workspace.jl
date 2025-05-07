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
    fast_ll_update!(a, groups_concerned)

    swap_node_labels!(a, s.u, s.v)
end


## below can be specialised for Bernoulli probably

function fast_ll_update!(a, groups_concerned)
    for g in groups_concerned
        a.log_likelihood[g[1], g[2]] = _fast_ll_one_group(a, g[1], g[2])
    end
end


function _fast_ll_one_group(a::Assignment, g1, g2)
    nodes_g1 = findall(x -> x == g1, a.node_labels)
    nodes_g2 = findall(x -> x == g2, a.node_labels)
    ll = 0.0
    d = a.θ[g1, g2]
    for u in nodes_g1
        for (v,e) in iterate_neighbors(a.edges,u) # assume implicitly that g1 != g2
            if v in nodes_g2
                if (g1 == g2 && u < v) || g1 != g2
                    ll += loglikelihood(d, e)
                end
            end
        end
    end
    return ll
end
