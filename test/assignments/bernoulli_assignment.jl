import NetworkHistogram as NH

@testset "test construction Bernoulli assignment" begin
    using Distributions: Bernoulli
    A = [0 1 1 1 0 0 1 0
         1 0 1 1 0 0 0 0
         1 1 0 0 0 0 0 0
         1 1 0 0 0 0 0 1
         0 0 0 0 0 1 1 1
         0 0 0 0 1 0 1 1
         1 0 0 0 1 1 0 0
         0 0 0 1 1 1 0 0]
    obs = NH.Observations(A, Bernoulli(0.5))
    node_labels = [1, 1, 1, 1, 2, 2, 2, 2]
    group_size = NH.GroupSize(8, 4)
    a = NH.BernoulliAssignment(obs, group_size, node_labels)
    for i in 1:8
        @test NH.get_group_of_vertex(a, i) == node_labels[i]
    end
    @test all(a.additional_data.A .== A)
    @test a.additional_data.realized == [5 2; 2 5]
    @test a.additional_data.counts == [6 16; 16 6]
    @test a.additional_data.estimated_theta == [5/6 1/8; 1/8 5/6]
end

@testset "test Bernoulli swap" begin
    using ..TestNetworkHistogram: test_swap_revertible
    using Distributions: Bernoulli
    A = [0 1 1 1 0 0 1 0
         1 0 1 1 0 0 0 0
         1 1 0 0 0 0 0 0
         1 1 0 0 0 0 0 1
         0 0 0 0 0 1 1 1
         0 0 0 0 1 0 1 1
         1 0 0 0 1 1 0 0
         0 0 0 1 1 1 0 0]
    obs = NH.Observations(A, Bernoulli(0.5))
    a = NH.BernoulliAssignment(
        obs, NH.GroupSize(8, 4), [1, 1, 1, 1, 2, 2, 2, 2])
    swap = NH.make_swap(a, (1, 2))
    test_swap_revertible(a, swap, obs)
end
