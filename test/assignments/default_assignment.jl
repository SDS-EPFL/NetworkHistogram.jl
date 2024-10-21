import NetworkHistogram as NH

@testset "test default swap" begin
    using ..TestNetworkHistogram: test_swap_revertible
    import Random, LinearAlgebra
    using Distributions: Bernoulli, Normal
    Random.seed!(1234)
    n = 20
    k = 5
    #data = LinearAlgebra.Symmetric(Random.rand(Bool,n,n))
    data = Random.rand(Normal(), n, n)
    g = NH.Observations(data, Normal(0, 1))
    labels = repeat(1:(n ÷ k), inner = k)
    a = NH.Assignment(NH.GroupSize(n, k), labels)
    swap = NH.DefaultSwap(1, 2)
    test_swap_revertible(a, swap, g)
end
