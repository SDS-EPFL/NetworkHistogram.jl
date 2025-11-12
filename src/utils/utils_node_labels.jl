function ordered_start_labels(n::Int, k::Int)
    labels = Vector{Int}(undef, n)
    base_size = n ÷ k
    remainder = n % k
    for group in 1:k
        fill!(view(labels, ((group - 1) * base_size + 1):(group * base_size)), group)
    end
    if remainder > 0
        fill!(view(labels, (k * base_size + 1):(k * base_size + remainder)), k)
    end
    return labels
end

function align_res_true_latents!(res::NethistResult, latents; type = :greedy)
    if type ∉ (:opt, :greedy)
        error("Unknown alignment type: $type. Use :opt or :greedy.")
    end
    if type == :opt
        @warn "The :opt alignment may not work that well; consider using :greedy instead."
    end
    new_labels, mapping = order_groups(res.labels, latents, Val(type))
    res.labels .= new_labels
    perm = [key for (key, val) in sort(collect(mapping), by = last)]
    permute!(res.model, perm)
end

function permute!(sbm, perm)
    permuted_theta = copy(sbm.θ)
    sbm.θ .= permuted_theta[perm, perm]
    sbm.size .= sbm.size[perm]
    sbm.cumsize .= cumsum(sbm.size)
end

function order_groups(node_labels, latents::AbstractVector, ::Val{:opt})
    return align_partitions(node_labels, latents)
end

function order_groups(node_labels, latents::AbstractVector, ::Val{:greedy})
    n = length(node_labels)
    k = length(unique(node_labels))
    sort_perm = sortperm(latents)
    sorted_group_labels = node_labels[sort_perm]
    dummy_group_labels = repeat(1:k, inner = n ÷ k + 1)[1:n]
    counts = Dict(group => countmap(dummy_group_labels[sorted_group_labels .== group])
    for group in 1:k)
    perm = sort(
        1:k, by = x -> Tuple(get(counts[x], g, 0) for g in 1:k), rev = true)
    new_labels = map(x -> findfirst(==(x), perm), node_labels)
    mapping = Dict(perm[i] => i for i in 1:k)
    return new_labels, mapping
end

function get_num_obs(A::AbstractMatrix)
    n = size(A, 1)
    return n * (n - 1) ÷ 2
end
