import NetworkHistogram as NH

@testset "test conversion to categorical observations" begin end

@testset "test Categorical swap" begin
    using ..TestNetworkHistogram: test_swap_revertible
    using Distributions: Categorical
    using LinearAlgebra: Symmetric
    import Random
    m = 5
    p = ones(m) ./ m
    n = 12
    k = 3
    dist = Categorical(p)
    sbm = NH.initialize_sbm(ones(k) ./ k, dist)
    A, _ = NH.sample(sbm, repeat(1:k, inner = n ÷ k))
    obs = NH.Observations(collect(A), dist)
    node_labels = repeat(1:k, inner = n ÷ k)
    a = NH.CategoricalAssignment(obs, NH.GroupSize(n, n ÷ k), node_labels)
    swap = NH.make_swap(a, (1, k + 1))
    test_swap_revertible(a, swap, obs)
end
