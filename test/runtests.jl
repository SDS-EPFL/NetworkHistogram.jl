using NetworkHistogram
using Test

using JLD
include("simple_test_example.jl")

@testset "NetworkHistogram.jl" begin
    include("pipeline_test.jl")
    include("test_multilayer.jl")
    include("proposal_test.jl")
    include("starting_labels_test.jl")
    include("oracle_bandwidth_test.jl")
    include("error_handling_tests.jl")
    include("config_rules/config_rule_test.jl")
    include("group_numbering_tests.jl")
end
