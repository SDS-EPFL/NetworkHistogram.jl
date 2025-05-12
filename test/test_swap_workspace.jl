using Test
using NetworkHistogram
using StatsBase
using Random


function manual_loglikelihood(A, node_labels, θ)
    n = size(A, 1)
    k = size(θ, 1)
    ll = 0.0
    for j in 1:n
        for i in 1:n
            if i!=j
                g1 = node_labels[i]
                g2 = node_labels[j]
                ll += NetworkHistogram.logpdf(θ[g1,g2], A[i,j])
            end
        end
    end
    return ll/2
end

function slow_swap(a::NetworkHistogram.Assignment, s::NetworkHistogram.Swap)
    labels = deepcopy(a.node_labels)
    labels[s.u], labels[s.v] = labels[s.v], labels[s.u]
    return NetworkHistogram.Assignment(labels, a.edges, a.θ[1,1])
end


@testset "Swap workspace likelihood update (Bernoulli)" begin
    Random.seed!(42)
    n = 6
    k = 2
    p1, p2 = 0.8, 0.3
    d = NetworkHistogram.Bernoulli(0.5)
    # Create a block model with two groups
    sbm = NetworkHistogram.BlockModel(k, d)
    sbm[1,1] = NetworkHistogram.Bernoulli(p1)
    sbm[2,2] = NetworkHistogram.Bernoulli(p2)
    sbm[1,2] = NetworkHistogram.Bernoulli(0.1)

    labels = StatsBase.inverse_rle(1:k, fill(n÷k, k))
    A = NetworkHistogram.sample(sbm, labels)
    edgelist = NetworkHistogram.EdgeList(A)
    assignment = NetworkHistogram.Assignment(labels, edgelist, NetworkHistogram.Dist(d))

    ll_original = NetworkHistogram.loglikelihood(assignment)
    ll_manual = manual_loglikelihood(A, assignment.node_labels, assignment.θ)
    @test isapprox(ll_original, ll_manual; atol=1e-10)

    # Swap two nodes from different groups
    indices = (1, n)
    swap = NetworkHistogram.make_swap(assignment, indices)
    slow_swapped = slow_swap(assignment, swap)
    NetworkHistogram.apply_swap!(assignment, swap)
    ll_after_swap =  NetworkHistogram.loglikelihood(assignment)
    ll_slow_swap = NetworkHistogram.loglikelihood(slow_swapped)
    ll_manual_after_swap = manual_loglikelihood(A, assignment.node_labels, assignment.θ)
    @test isapprox(ll_after_swap, ll_manual_after_swap; atol=1e-10)
    @test isapprox(ll_after_swap, ll_slow_swap; atol=1e-10)

    # Revert the swap
    NetworkHistogram.revert_swap!(assignment, swap)
    ll_after_revert = NetworkHistogram.loglikelihood(assignment)
    ll_manual_after_revert = manual_loglikelihood(A, assignment.node_labels, assignment.θ)
    @test isapprox(ll_after_revert, ll_manual_after_revert; atol=1e-10)
    @test isapprox(ll_after_revert, ll_original; atol=1e-10)
end
