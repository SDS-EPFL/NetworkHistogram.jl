using Test
using LinearAlgebra, SparseArrays
using NetworkHistogram
@testset "Tests" begin

    include("test_data_format.jl")
    include("test_distributions_type.jl")
    include("test_swap_workspace.jl")
end
