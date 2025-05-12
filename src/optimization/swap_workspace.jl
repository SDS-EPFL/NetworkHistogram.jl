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


function apply_swap!(a::Assignment, s::Swap)
    # swap node labels
    swap_node_labels!(a, s.u, s.v)
    new_assignment = Assignment(a.node_labels, a.edges, a.θ[1,1])
    a.θ = new_assignment.θ
    a.log_likelihood = new_assignment.log_likelihood
    # # fully rebuild θ and log_likelihood based on new labels
    # k = size(a.θ, 1)
    # # initial distribution template and zero-likelihood
    # base_dist = a.θ[1, 1]
    # a.θ = SymArray(k, base_dist)
    # a.log_likelihood = SymArray(k, zero(eltype(a.log_likelihood)))
    # # accumulate edge contributions
    # for u in 1:length(a.node_labels)
    #     g_u = group(a, u)
    #     for (v, d) in iterate_neighbors(a.dists, u)
    #         if u < v
    #             g_v = group(a, v)
    #             a.θ[g_u, g_v] = add_to(a.θ[g_u, g_v], d)
    #         end
    #     end
    # end
    # # recompute log likelihoods for all group pairs
    # for g1 in 1:k, g2 in g1:k
    #     edges = get_edges_in_groups(a, g1, g2)
    #     a.log_likelihood[g1, g2] = loglikelihood(a.θ[g1, g2], edges)
    # end
end


## below can be specialised for Bernoulli probably

function fast_ll_update!(a, groups_concerned)
    for g in groups_concerned
        # Use get_edges_in_groups to get the correct set of edges
        edges = get_edges_in_groups(a, g[1], g[2])
        a.log_likelihood[g[1], g[2]] = loglikelihood(a.θ[g[1], g[2]], edges)
    end
end
