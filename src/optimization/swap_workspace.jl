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

apply_swap!(a::Assignment, s::Swap) = _slow_swap!(a, s)


## below is not faster than the above, need to find a way to take advantage of the sparsity
# somewhere in the datastructure for it to make a difference

# function apply_swap!(a::Assignment, s::Swap)
#     g1 = group(a, s.u)
#     g2 = group(a, s.v)
#     groups_concerned = Set([minmax(g1, g2)])
#     for (u, g_old, g_new) in [(s.u, g1, g2), (s.v, g2, g1)]
#         for (v,d) in iterate_neighbors(a.dists, u)
#             if v == s.u || v == s.v
#                 continue
#             end
#             g_v = group(a, v)
#             a.θ[g_new, g_v] = add_to(a.θ[g_new, g_v], d)
#             a.θ[g_old, g_v] = remove_from(a.θ[g_old, g_v], d)
#             push!(groups_concerned, minmax(g_new, g_v), minmax(g_old, g_v))
#         end
#     end
#     swap_node_labels!(a, s.u, s.v)
#     fast_ll_update!(a, groups_concerned)
# end


# ## below can be specialised for Bernoulli probably (probably above needs to be actually)

# function fast_ll_update!(a, groups_concerned)
#     for g in groups_concerned
#         a.log_likelihood[g[1], g[2]] = _fast_ll_one_group(a, g[1], g[2])
#     end
# end


# function _fast_ll_one_group(a::Assignment, g1, g2)
#     nodes_g1 = findall(x -> x == g1, a.node_labels)
#     nodes_g2 = findall(x -> x == g2, a.node_labels)
#     ll = 0.0
#     d = a.θ[g1, g2]
#     for u in nodes_g1
#         for (v,e) in iterate_neighbors(a.edges,u) # assume implicitly that g1 != g2
#             if v in nodes_g2
#                 if (g1 == g2 && u < v) || g1 != g2
#                     ll += logpdf(d, e)
#                 end
#             end
#         end
#     end
#     return ll
# end
