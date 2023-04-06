using NetworkHistogram
using Test, RCall

using JLD

@testset "NetworkHistogram.jl" begin
    include("pipeline_test.jl")
    include("proposal_test.jl")
    include("starting_labels_test.jl")
    include("oracle_bandwidth_test.jl")
end
