using Test
using LinearAlgebra, SparseArrays
using NetworkHistogram

@testset "Tests" begin
    include("test_data_format.jl")
    include("test_distributions_type.jl")
    include("test_swap_workspace.jl")
    include("test_cat_case.jl")
    include("test_get_edges_in_groups.jl")
end
