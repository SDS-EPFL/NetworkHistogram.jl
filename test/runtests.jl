using Test

include("TestNetworkHistogram.jl")

@testset "Assignment tests" begin
    include("assignments/default_assignment.jl")
    include("assignments/bernoulli_assignment.jl")
    include("assignments/categorical_assignment.jl")
end

@testset "Rule optimization tests" begin
    include("optimisation/config_rules/init_rule.jl")
end
