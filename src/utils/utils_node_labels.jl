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
    return align_partitions(node_labels, latents)[2]
end

function get_num_obs(A::AbstractMatrix)
    n = size(A, 1)
    return n * (n - 1) ÷ 2
end
