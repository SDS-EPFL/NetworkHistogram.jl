using NetworkHistogram
using Test

using JLD

@testset "NetworkHistogram.jl" begin
    include("pipeline_test.jl")
    include("proposal_test.jl")
    include("starting_labels_test.jl")
    include("error_handling_tests.jl")

    include("config_rules/config_rule_test.jl")
end
