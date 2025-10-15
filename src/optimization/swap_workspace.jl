mutable struct WorkspaceSwap{D, F}
    θ::SymArray{D}
    log_likelihood_per_group::SymArray{F}
end

function make_workspace(a::Assignment)
    # Pre-allocate workspace with same structure
    k = number_groups(a)
    θ_copy = SymArray(k, zero(a.θ[1, 1]))
    ll_copy = SymArray(k, 0.0)
    return WorkspaceSwap(θ_copy, ll_copy)
end

mutable struct Swap{W}
    u::Int
    v::Int
    workspace::W
end

function copy_symarray!(dest::SymArray, src::SymArray)
    # In-place copy without allocation
    # SymArray stores data in a dictionary .d
    # Just overwrite the values - don't empty first!
    @inbounds for key in keys(src.d)
        dest.d[key] = src.d[key]
    end
end

function make_swap_workspace!(ws, a::Assignment)
    # Use in-place copy instead of deepcopy
    copy_symarray!(ws.θ, a.θ)
    copy_symarray!(ws.log_likelihood_per_group, a.log_likelihood)
end

function revert_swap_workspace!(a::Assignment, ws)
    # Use in-place copy instead of deepcopy
    copy_symarray!(a.θ, ws.θ)
    copy_symarray!(a.log_likelihood, ws.log_likelihood_per_group)
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

# for reference and testing
function _slow_swap!(a::Assignment, s::Swap)
    swap_node_labels!(a, s.u, s.v)
    a.θ,
    a.log_likelihood = _compute_theta_and_ll(
        a.node_labels, a.dists, a.edges, a.θ[1, 1])
end

# apply_swap!(a::Assignment, s::Swap) = _slow_swap!(a, s)

function apply_swap!(a::Assignment, s::Swap)
    u, v = s.u, s.v
    gu = a.node_labels[u]
    gv = a.node_labels[v]

    # Pre-allocate with reasonable capacity to avoid resizing
    # Most swaps affect at most degree(u) + degree(v) + 1 group pairs
    groups_concerned = Set{Tuple{Int, Int}}()
    sizehint!(groups_concerned, 16)  # Reasonable default
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
