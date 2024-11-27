using Test
using Aqua

include("TestNetworkHistogram.jl")

@testset "Tests" begin
    @testset "Assignment tests" begin
        include("assignments/default_assignment.jl")
        include("assignments/bernoulli_assignment.jl")
        include("assignments/categorical_assignment.jl")
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
