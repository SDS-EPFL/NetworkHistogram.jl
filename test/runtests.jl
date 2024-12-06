using Test
using Aqua
using SparseArrays
include("TestNetworkHistogram.jl")

@testset "Tests" begin

    @testset "Discretizer tests" begin
        include("discretised_dist/discretizer.jl")
    end
    @testset "Assignment tests" begin
        include("assignments/default_assignment.jl")
        include("assignments/bernoulli_assignment.jl")
        include("assignments/categorical_assignment.jl")
        include("assignments/sparse_assignment.jl")
    end

    @testset "Rule optimization tests" begin
        include("optimisation/config_rules/init_rule.jl")
    end

    @testset "Observations tests" begin
        include("observations/discretisation.jl")
    end

    @testset "API tests" begin
        include("test_api.jl")
    end
    @testset "Aqua.jl for package quality" begin
        using NetworkHistogram
        Aqua.test_all(NetworkHistogram)
    end
end
