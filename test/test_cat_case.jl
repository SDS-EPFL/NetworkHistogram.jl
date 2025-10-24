using Test
using NetworkHistogram
using StatsBase
using Random
using Distributions
using StaticArrays

@testset "Swap workspace likelihood update (Categorical)" begin
    Random.seed!(42)
    n = 10
    k = 2
    m = 3
    ps = SVector{m}(fill(1 / m, m))
    d_mine = NetworkHistogram.Cat(ps)

    θ = [NetworkHistogram.Cat(SVector{3}([0.7, 0.2, 0.1])) NetworkHistogram.Cat(SVector{3}([0.1, 0.3, 0.6]));
         NetworkHistogram.Cat(SVector{3}([0.1, 0.3, 0.6])) NetworkHistogram.Cat(SVector{3}([0.3, 0.4, 0.3]))]
    sbm = DecoratedSBM(θ, [0.5, 0.5])

    labels = StatsBase.inverse_rle(1:k, fill(n ÷ k, k))
    latents = vcat(repeat([0.2], n ÷ 2), repeat([0.8], n ÷ 2))
    A = sample_graph(sbm, latents)

    edgelist = NetworkHistogram.EdgeList(A)
    assignment = NetworkHistogram.Assignment(
        labels, edgelist, NetworkHistogram.Dist(d_mine))

    for ind in eachindex(assignment.additional_workspace.counts)
        @test assignment.additional_workspace.counts[ind] ==
              sum(assignment.additional_workspace.realized[ind])
    end

    ll_original = NetworkHistogram.loglikelihood(assignment)

    # Swap two nodes from different groups
    indices = (1, n)
    swap = NetworkHistogram.make_swap(assignment, indices)
    true_swapped = deepcopy(labels)
    true_swapped[1] = labels[n]
    true_swapped[n] = labels[1]
    NetworkHistogram.apply_swap!(assignment, swap)

    nodes_label_swapped = deepcopy(assignment.node_labels)
    new_a = NetworkHistogram.Assignment(
        nodes_label_swapped, edgelist, NetworkHistogram.Dist(d_mine))
    ll_new_a = NetworkHistogram.loglikelihood(new_a)
    ll_after_swap = NetworkHistogram.loglikelihood(assignment)
    ws_new = new_a.additional_workspace
    ws_old = assignment.additional_workspace
    @test ws_new.counts == ws_old.counts
    @test ws_new.realized == ws_old.realized
    @test ws_new.estimated == ws_old.estimated

    @test new_a.node_labels == assignment.node_labels
    @test new_a.node_labels == true_swapped
    @test new_a.log_likelihood == assignment.log_likelihood
    @test isapprox(ll_after_swap, ll_new_a; atol = 1e-10)

    # Revert the swap
    NetworkHistogram.revert_swap!(assignment, swap)
    ll_after_revert = NetworkHistogram.loglikelihood(assignment)
    @test isapprox(ll_after_revert, ll_original; atol = 1e-10)
end
