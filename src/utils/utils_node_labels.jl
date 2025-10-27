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

function node_labels_to_latents(node_labels::AbstractVector{Int}, sbm)
    return map(label -> _label_to_latent(label, sbm), node_labels)
end

function _label_to_latent(label::Int, sbm)
    return sbm.cumsize[label] - eps()
end

function align_res_true_latents!(res::NethistResult, latents)
    perm = order_groups(res.labels, latents)
    permute!(res.model, perm)
    res.labels .= map(x -> findfirst(==(x), perm), res.labels)
end

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
    counts = Dict(group => countmap(dummy_group_labels[sorted_group_labels .== group])
    for group in 1:k)
    return sort(
        1:k, by = x -> Tuple(get(counts[x], g, 0) for g in 1:k), rev = true)
end
