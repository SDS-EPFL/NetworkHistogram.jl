import NetworkHistogram as NH

@testset "regression test" begin
    using Distributions: Bernoulli
    A = BitMatrix([0 0 1 0 1 0 1 1 0 1
                   0 0 1 1 1 1 1 1 0 0
                   1 1 0 1 0 0 0 0 1 0
                   0 1 1 0 1 0 1 0 0 0
                   1 1 0 1 0 0 1 0 0 1
                   0 1 0 0 0 0 0 1 0 0
                   1 1 0 1 1 0 0 1 0 1
                   1 1 0 0 0 1 1 0 0 1
                   0 0 1 0 0 0 0 0 0 1
                   1 0 0 0 1 0 1 1 1 0])
    h_true_nethist = 2.643731 # version 0.2.3 from nethist package
    k_true = 3
    obs = NH.Observations(A, Bernoulli(0.5))
    @testset "degrees" begin
        k = NH.select_number_node_per_block(obs, NH.EstimatedDegrees())
        @test k == k_true
    end
    @testset "eigenvalues" begin
        k = NH.select_number_node_per_block(obs, NH.EstimatedEigenvalues())
        @test k == k_true
    end
end

@testset "test oracle K" begin
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
    oracle = NH.OracleK(4)
    @test NH.select_number_node_per_block(obs, oracle) == 4
    err = ArgumentError("The number of blocks 5 is too large for the number of nodes \
        8, it should be at most 4")
    @test_throws err NH.select_number_node_per_block(obs, NH.OracleK(5))
end
