import NetworkHistogram as NH

using Random


@testset "test sparse swap" begin
    Random.seed!(1234123)
    using ..TestNetworkHistogram: test_swap_revertible, to_default_assignment
    using Distributions: Categorical
    using LinearAlgebra: Symmetric
    import Random
    m = 4
    p = ones(m) ./ m
    n = 12
    k = 4
    dist = Categorical(p)
    sbm = NH.initialize_sbm(ones(k) ./ k, dist)
    node_labels = repeat(1:k, inner = n ÷ k)
    A = sparse(first(NH.sample(sbm, node_labels)))
    g = NH.Observations(A, dist)
    a = NH.SparseAssignment(g, NH.GroupSize(n, n ÷ k), node_labels)
    swap = NH.make_swap(a, (1, k + 1))
    @test A[:, 1] != A[:, k + 1]
    a_test = deepcopy(a)
    NH.apply_swap!(a_test, swap)
    @test NH.get_group_of_vertex(a, swap.index1) ==
          NH.get_group_of_vertex(a_test, swap.index2)
    @test NH.get_group_of_vertex(a, swap.index2) ==
          NH.get_group_of_vertex(a_test, swap.index1)
    # force recomputation of the log likelihood using default assignment
    a_new = to_default_assignment(a_test)
    @test NH.loglikelihood(a_new, g) ≈ NH.loglikelihood(a_test, g)
    @test a_test.additional_data.realized != a.additional_data.realized
    @test a_test.additional_data.estimated_theta !=
          a.additional_data.estimated_theta
    @test a_test.additional_data.log_likelihood !=
          a.additional_data.log_likelihood
    # revert the swap and check if the assignment is the same as before
    NH.revert_swap!(a_test, swap)
    @test a == a_test
    @test NH.loglikelihood(a, g) ≈ NH.loglikelihood(a_test, g)
end

@testset "fast sparse update test" begin
    using Distributions
    realized = [[[1, 0, 0]] [[0, 4, 0]] [[0, 0, 4]];
                [[0, 4, 0]] [[1, 0, 0]] [[0, 0, 4]];
                [[0, 0, 4]] [[0, 0, 4]] [[1, 0, 0]]]
    realized = [realized[I][k]
                for k in eachindex(realized[1, 1]),
    I in CartesianIndices(realized)]
    counts = [1 4 4
              4 1 4
              4 4 1]
    A = sparse([0 1 2 2 3 3
         1 0 2 2 3 3
         2 2 0 1 3 3
         2 2 1 0 3 3
         3 3 3 3 0 1
         3 3 3 3 1 0])
    groupsize = NH.GroupSize(6, 2)
    node_labels = [1, 1, 2, 2, 3, 3]
    g = NH.Observations(A, Categorical(3))
    k = 3
    m = 3
    n = size(A,1)
    a = NH.SparseAssignment(g, NH.GroupSize(n, n ÷ k), node_labels)
    for index in eachindex(realized)
        @test all(realized[index] .== a.additional_data.realized[index])
    end
    @test loglikelihood(a, g) ≈ 0
    @test a.additional_data.counts == counts
    swap_id = (1, 3)
    ras = [[[0, 1, 0]] [[2, 2, 0]] [[0, 0, 4]];
           [[2, 2, 0]] [[0, 1, 0]] [[0, 0, 4]];
           [[0, 0, 4]] [[0, 0, 4]] [[1, 0, 0]]]
    realized_after_swap = [ras[I][k]
                           for k in eachindex(ras[1, 1]),
    I in CartesianIndices(ras)]

    swap = NH.make_swap(a, swap_id)
    NH.apply_swap!(a, swap)
    for j in 1:3
        for i in 1:3
            @test all(realized_after_swap[:, i, j] .==
                      a.additional_data.realized[:, i, j])
            @test all(a.additional_data.estimated_theta[:, i, j] .≈
                      realized_after_swap[:, i, j] ./ counts[i, j])
        end
    end
    @test loglikelihood(a, g) == 4 * log(0.5)
end
