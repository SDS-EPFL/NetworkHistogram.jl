# Reasonable default capacity for affected groups in a swap
const MAX_AFFECTED_GROUPS = 16

mutable struct WorkspaceSwap{D, F, G}
    θ::SymArray{D}
    log_likelihood_per_group::SymArray{F}
    groups_buffer::G  # Pre-allocated buffer for affected group pairs
end

function make_workspace(a::Assignment)
    # Pre-allocate workspace with same structure
    k = number_groups(a)
    θ_copy = SymArray(k, zero(a.θ[1, 1]))
    ll_copy = SymArray(k, 0.0)
    groups_buffer = Set{Tuple{Int, Int}}()
    sizehint!(groups_buffer, MAX_AFFECTED_GROUPS)
    return WorkspaceSwap(θ_copy, ll_copy, groups_buffer)
end

mutable struct Swap{W}
    u::Int
    v::Int
    workspace::W
end

function make_swap_workspace!(ws, a::Assignment)
    # Use in-place copy instead of deepcopy
    copy!(ws.θ, a.θ)
    copy!(ws.log_likelihood_per_group, a.log_likelihood)
end

function revert_swap_workspace!(a::Assignment, ws)
    # Use in-place copy instead of deepcopy
    copy!(a.θ, ws.θ)
    copy!(a.log_likelihood, ws.log_likelihood_per_group)
end

function make_swap(a::Assignment, id)
    ws = make_workspace(a)
    make_swap_workspace!(ws, a)  # Actually copy the current state
    return Swap(id[1], id[2], ws)
end

function make_swap!(swap::Swap, a::Assignment, id)
    swap.u, swap.v = id
    make_swap_workspace!(swap.workspace, a)
end

function revert_swap!(assignment::Assignment, swap::Swap)
    # swap labels back to original
    swap_node_labels!(assignment, swap.u, swap.v)
    # restore saved θ and log likelihoods
    revert_swap_workspace!(assignment, swap.workspace)
end

function swap_node_labels!(a::Assignment, i, j)
    a.node_labels[i], a.node_labels[j] = a.node_labels[j], a.node_labels[i]
end

function apply_swap!(a::Assignment, s::Swap)
    u, v = s.u, s.v
    gu = a.node_labels[u]
    gv = a.node_labels[v]

    # Reuse pre-allocated buffer instead of allocating new Set each time
    groups_concerned = s.workspace.groups_buffer
    empty!(groups_concerned)
    push!(groups_concerned, minmax(gu, gv))

    @inbounds for (node, d) in iterate_neighbors(a.dists, u)
        if node == v
            continue
        end
        g1 = a.node_labels[node]
        a.θ[gv, g1] = add_to(a.θ[gv, g1], d)
        a.θ[gu, g1] = remove_from(a.θ[gu, g1], d)
        push!(groups_concerned, minmax(gu, g1))
        push!(groups_concerned, minmax(gv, g1))
    end

    @inbounds for (index, (node, d)) in enumerate(iterate_neighbors(a.dists, v))
        if node == u
            continue
        end
        g2 = a.node_labels[node]
        a.θ[gu, g2] = add_to(a.θ[gu, g2], d)
        a.θ[gv, g2] = remove_from(a.θ[gv, g2], d)
        push!(groups_concerned, minmax(gv, g2))
        push!(groups_concerned, minmax(gu, g2))
    end

    swap_node_labels!(a, u, v)
    @inbounds for (g1, g2) in groups_concerned
        a.log_likelihood[g1, g2] = 0.0
        for e in get_edges_in_groups(a.node_labels, a.edges, g1, g2)
            a.log_likelihood[g1, g2] += logpdf(a.θ[g1, g2], e)
        end
    end
end
