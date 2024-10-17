using Test

include("TestNetworkHistogram.jl")

@testset "Assignment tests" begin
    include("assignments/default_assignment.jl")
    include("assignments/bernoulli_assignment.jl")
    include("assignments/categorical_assignment.jl")
end
