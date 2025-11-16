function ordered_start_labels(n::Int, k::Int)
    labels = Vector{Int}(undef, n)
    base_size = n ÷ k
    remainder = n % k
    for group = 1:k
        fill!(view(labels, ((group-1)*base_size+1):(group*base_size)), group)
    end
    if remainder > 0
        fill!(view(labels, (k*base_size+1):(k*base_size+remainder)), k)
    end
    return labels
end

function align_res_true_latents!(res::NethistResult, latents)
    new_labels, mapping = order_groups(res.labels, latents)
    res.labels .= new_labels
    perm = [key for (key, val) in sort(collect(mapping), by = last)]
    permute!(res.model, perm)
end

#TODO: move to Graphons.jl see https://github.com/SDS-EPFL/Graphons.jl/pull/17
function permute!(sbm, perm)
    permuted_theta = copy(sbm.θ)
    sbm.θ .= permuted_theta[perm, perm]
    sbm.size .= sbm.size[perm]
    sbm.cumsize .= cumsum(sbm.size)
end

function order_groups(node_labels, latents::AbstractVector)
    n = length(node_labels)
    k = length(unique(node_labels))
    sort_perm = sortperm(latents)
    sorted_group_labels = node_labels[sort_perm]
    dummy_group_labels = repeat(1:k, inner = n ÷ k + 1)[1:n]
    counts = Dict(
        group => countmap(dummy_group_labels[sorted_group_labels .== group]) for
        group = 1:k
    )
    perm = sort(1:k, by = x -> Tuple(get(counts[x], g, 0) for g = 1:k), rev = true)
    new_labels = map(x -> findfirst(==(x), perm), node_labels)
    mapping = Dict(perm[i] => i for i = 1:k)
    return new_labels, mapping
end

function get_num_obs(A::AbstractMatrix)
    n = size(A, 1)
    return n * (n - 1) ÷ 2
end

"""
Align the source and target matrices using optimal transport. This function requires
the PythonCall.jl package to be loaded
"""
function align_matrices end

"""
Get the permutation aligning source and target matrices using optimal transport. This function requires
the PythonCall.jl package to be loaded
"""
function get_perm_alignment end
