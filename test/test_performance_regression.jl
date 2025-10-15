using Test
using NetworkHistogram
using StatsBase
using Random
using Distributions
using StaticArrays
using BenchmarkTools

"""
Performance regression test suite for NetworkHistogram optimization.

This file contains benchmarks for the key optimization operations, designed to:
1. Track performance improvements/regressions over time
2. Identify bottlenecks in the optimization workflow
3. Ensure optimization changes maintain correctness

Based on the workflow in test_decorated_paper.jl
"""

# Helper function to create test networks
function create_test_sbm_bernoulli(n_groups::Int, n_nodes::Int; seed = 42)
    Random.seed!(seed)
    d = NetworkHistogram.Bernoulli(0.5)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    # Create varied probabilities between groups
    for g1 in 1:n_groups
        for g2 in g1:n_groups
            p = 0.1 + 0.7 * rand()
            sbm[g1, g2] = NetworkHistogram.Bernoulli(p)
        end
    end

    # Ensure labels has exactly n_nodes elements
    base_size = n_nodes ÷ n_groups
    remainder = n_nodes % n_groups
    sizes = fill(base_size, n_groups)
    sizes[1:remainder] .+= 1  # Distribute remainder to first groups
    labels = StatsBase.inverse_rle(1:n_groups, sizes)
    @assert length(labels) == n_nodes
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

function create_test_sbm_categorical(
        n_groups::Int, n_nodes::Int, n_categories::Int; seed = 42)
    Random.seed!(seed)
    ps = SVector{n_categories}(fill(1 / n_categories, n_categories))
    d = NetworkHistogram.Cat(ps)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    # Create varied probability distributions between groups
    for g1 in 1:n_groups
        for g2 in g1:n_groups
            probs = rand(n_categories)
            probs ./= sum(probs)
            sbm[g1, g2] = NetworkHistogram.Cat(SVector{n_categories}(probs))
        end
    end

    # Ensure labels has exactly n_nodes elements
    base_size = n_nodes ÷ n_groups
    remainder = n_nodes % n_groups
    sizes = fill(base_size, n_groups)
    sizes[1:remainder] .+= 1  # Distribute remainder to first groups
    labels = StatsBase.inverse_rle(1:n_groups, sizes)
    @assert length(labels) == n_nodes
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

@testset "Performance Regression Tests" begin
    @testset "Bernoulli Networks" begin
        @testset "Small network (n=50, k=2)" begin
            A, labels, d = create_test_sbm_bernoulli(2, 50)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            # Benchmark single swap operation
            swap = NetworkHistogram.make_swap(assignment, (1, 50))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end setup=(NetworkHistogram.make_swap_workspace!($swap.workspace, $assignment)) samples=100 evals=1

            # Verify correctness
            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Bernoulli (n=50, k=2) - Single swap" median=median(b_swap.times) / 1e6 mean=mean(b_swap.times) /
                                                                                               1e6
        end

        @testset "Medium network (n=200, k=3)" begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            swap = NetworkHistogram.make_swap(assignment, (1, 200))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end setup=(NetworkHistogram.make_swap_workspace!($swap.workspace, $assignment)) samples=100 evals=1

            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Bernoulli (n=200, k=3) - Single swap" median=median(b_swap.times) / 1e6 mean=mean(b_swap.times) /
                                                                                                1e6
        end

        @testset "Large network (n=500, k=5)" begin
            A, labels, d = create_test_sbm_bernoulli(5, 500)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            swap = NetworkHistogram.make_swap(assignment, (1, 500))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end setup=(NetworkHistogram.make_swap_workspace!($swap.workspace, $assignment)) samples=50 evals=1

            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Bernoulli (n=500, k=5) - Single swap" median=median(b_swap.times) / 1e6 mean=mean(b_swap.times) /
                                                                                                1e6
        end
    end

    @testset "Categorical Networks" begin
        @testset "Small network (n=50, k=2, m=3)" begin
            A, labels, d = create_test_sbm_categorical(2, 50, 3)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            swap = NetworkHistogram.make_swap(assignment, (1, 50))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end samples=100 evals=1

            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Categorical (n=50, k=2, m=3) - Single swap" median=median(b_swap.times) /
                                                                      1e6 mean=mean(b_swap.times) /
                                                                               1e6
        end

        @testset "Medium network (n=200, k=3, m=4)" begin
            A, labels, d = create_test_sbm_categorical(3, 200, 4)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            swap = NetworkHistogram.make_swap(assignment, (1, 200))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end samples=100 evals=1

            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Categorical (n=200, k=3, m=4) - Single swap" median=median(b_swap.times) /
                                                                       1e6 mean=mean(b_swap.times) /
                                                                                1e6
        end

        @testset "Large network (n=500, k=5, m=5)" begin
            A, labels, d = create_test_sbm_categorical(5, 500, 5)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            swap = NetworkHistogram.make_swap(assignment, (1, 500))
            ll_before = NetworkHistogram.loglikelihood(assignment)

            b_swap = @benchmark begin
                NetworkHistogram.apply_swap!($assignment, $swap)
                NetworkHistogram.revert_swap!($assignment, $swap)
            end samples=50 evals=1

            ll_after = NetworkHistogram.loglikelihood(assignment)
            @test isapprox(ll_before, ll_after, atol = 1e-10)

            @info "Categorical (n=500, k=5, m=5) - Single swap" median=median(b_swap.times) /
                                                                       1e6 mean=mean(b_swap.times) /
                                                                                1e6
        end
    end

    @testset "Full Optimization Workflow" begin
        @testset "Bernoulli - Short optimization (n=100, k=3)" begin
            A, labels, d = create_test_sbm_bernoulli(3, 100)

            # Randomize initial labels
            initial_labels = rand(1:3, 100)

            b_optimize = @benchmark begin
                # Create fresh params for each benchmark iteration
                params = NetworkHistogram.GreedyParams(
                    1_000,  # Small number for testing
                    NetworkHistogram.RandomNodeSwap(),
                    NetworkHistogram.Strict(),
                    NetworkHistogram.PreviousBestValue(500),
                    false  # No progress bar for benchmarking
                )
                NetworkHistogram.nethist($A, $d, $initial_labels, params)
            end samples=10 evals=1

            @info "Bernoulli full optimization (n=100, 1k iters)" median=median(b_optimize.times) /
                                                                         1e6 mean=mean(b_optimize.times) /
                                                                                  1e6
        end

        @testset "Categorical - Short optimization (n=100, k=3, m=3)" begin
            A, labels, d = create_test_sbm_categorical(3, 100, 3)

            # Randomize initial labels
            initial_labels = rand(1:3, 100)

            b_optimize = @benchmark begin
                # Create fresh params for each benchmark iteration
                params = NetworkHistogram.GreedyParams(
                    1_000,
                    NetworkHistogram.RandomNodeSwap(),
                    NetworkHistogram.Strict(),
                    NetworkHistogram.PreviousBestValue(500),
                    false
                )
                NetworkHistogram.nethist($A, $d, $initial_labels, params)
            end samples=10 evals=1

            @info "Categorical full optimization (n=100, 1k iters)" median=median(b_optimize.times) /
                                                                           1e6 mean=mean(b_optimize.times) /
                                                                                    1e6
        end
    end

    @testset "Component Benchmarks" begin
        @testset "Assignment creation (n=200, k=3)" begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)

            b_assignment = @benchmark begin
                NetworkHistogram.Assignment($labels, $edgelist, NetworkHistogram.Dist($d))
            end samples=100

            @info "Assignment creation (n=200)" median=median(b_assignment.times) / 1e6 mean=mean(b_assignment.times) /
                                                                                             1e6
        end

        @testset "EdgeList creation (n=200)" begin
            A, _, _ = create_test_sbm_bernoulli(3, 200)

            b_edgelist = @benchmark begin
                NetworkHistogram.EdgeList($A)
            end samples=100

            @info "EdgeList creation (n=200)" median=median(b_edgelist.times) / 1e6 mean=mean(b_edgelist.times) /
                                                                                         1e6
        end

        @testset "Loglikelihood computation (n=200, k=3)" begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            b_ll = @benchmark begin
                NetworkHistogram.loglikelihood($assignment)
            end samples=1000

            @info "Loglikelihood computation (n=200)" median=median(b_ll.times) / 1e3 mean=mean(b_ll.times) /
                                                                                           1e3
        end

        @testset "Get edges in groups (n=200, k=3)" begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))

            b_get_edges = @benchmark begin
                NetworkHistogram.get_edges_in_groups($assignment, 1, 2)
            end samples=1000

            @info "Get edges in groups (n=200)" median=median(b_get_edges.times) / 1e3 mean=mean(b_get_edges.times) /
                                                                                            1e3
        end
    end
end
