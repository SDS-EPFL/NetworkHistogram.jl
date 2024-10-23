import NetworkHistogram as NH

@testset "test conversion to categorical observations" begin end

@testset "test Categorical swap" begin
    using ..TestNetworkHistogram: test_swap_revertible, to_default_assignment
    using Distributions: Categorical
    using LinearAlgebra: Symmetric
    import Random
    m = 2
    p = ones(m) ./ m
    n = 12
    k = 4
    dist = Categorical(p)
    sbm = NH.initialize_sbm(ones(k) ./ k, dist)
    A, _ = NH.sample(sbm, repeat(1:k, inner = n ÷ k))
    g = NH.Observations(collect(A), dist)
    node_labels = repeat(1:k, inner = n ÷ k)
    a = NH.CategoricalAssignment(g, NH.GroupSize(n, n ÷ k), node_labels)
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

@testset "fast update test" begin
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
    A = [0 1 2 2 3 3
         1 0 2 2 3 3
         2 2 0 1 3 3
         2 2 1 0 3 3
         3 3 3 3 0 1
         3 3 3 3 1 0]
    groupsize = NH.GroupSize(6, 2)
    node_labels = [1, 1, 2, 2, 3, 3]
    g = NH.Observations(A, Categorical(3))
    a = NH.CategoricalAssignment(g, groupsize, node_labels)
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

#todo: test ll against categorical likelihood on basic assignment
@testset "test swap is not overwritten" begin
    A = [0 4 4 2 1 2 2 3 4 2 3 1 4 1 1 3 4 4 3 3
         4 0 4 2 4 2 1 1 1 3 3 1 1 1 3 3 4 2 1 4
         4 4 0 1 2 4 2 2 1 3 2 3 1 2 3 2 3 4 1 1
         2 2 1 0 2 1 2 2 2 3 1 1 3 3 3 3 3 1 1 2
         1 4 2 2 0 4 1 4 3 2 4 3 4 3 1 3 1 1 1 3
         2 2 4 1 4 0 2 3 1 3 1 4 3 3 1 3 1 3 3 3
         2 1 2 2 1 2 0 3 2 2 1 1 1 3 3 1 1 3 1 1
         3 1 2 2 4 3 3 0 4 3 2 3 1 1 1 1 1 3 2 1
         4 1 1 2 3 1 2 4 0 3 1 1 1 3 2 1 3 1 4 1
         2 3 3 3 2 3 2 3 3 0 1 3 1 1 3 1 3 1 1 4
         3 3 2 1 4 1 1 2 1 1 0 2 3 2 2 1 2 2 1 3
         1 1 3 1 3 4 1 3 1 3 2 0 4 4 2 2 2 3 1 1
         4 1 1 3 4 3 1 1 1 1 3 4 0 2 2 1 2 1 1 3
         1 1 2 3 3 3 3 1 3 1 2 4 2 0 1 2 1 2 1 1
         1 3 3 3 1 1 3 1 2 3 2 2 2 1 0 2 1 2 1 1
         3 3 2 3 3 3 1 1 1 1 1 2 1 2 2 0 1 1 1 3
         4 4 3 3 1 1 1 1 3 3 2 2 2 1 1 1 0 1 1 1
         4 2 4 1 1 3 3 3 1 1 2 3 1 2 2 1 1 0 1 1
         3 1 1 1 1 3 1 2 4 1 1 1 1 1 1 1 1 1 0 1
         3 4 1 2 3 3 1 1 1 4 3 1 3 1 1 3 1 1 1 0]
    g = NH.Observations(A, Categorical(4))
    h = 6
    a = NH.make_assignment(
        g, h, NH.InitRule(NH.OrderedStart(), Val{NH.CategoricalData}()))
    a_ref = deepcopy(a)
    swap_indices = [(18, 5), (15, 10), (5, 13)]
    swap = NH.make_swap(a, swap_indices[1])
    for swap_index in swap_indices
        NH.make_swap!(swap, a, swap_index)
        NH.apply_swap!(a, swap)
        @test swap.realized == a_ref.additional_data.realized
        @test swap.estimated_theta == a_ref.additional_data.estimated_theta
        NH.revert_swap!(a, swap)
    end
end
