using Test
using LinearAlgebra, SparseArrays
using NetworkHistogram

# Check if BenchmarkTools is available (it's not required for basic tests)
const RUN_BENCHMARKS = try
    using BenchmarkTools
    true
catch
    @warn "BenchmarkTools not available, skipping performance regression tests"
    false
end

@testset "Tests" begin
    include("test_data_format.jl")
    include("test_distributions_type.jl")
    include("test_swap_workspace.jl")
    include("test_cat_case.jl")
    include("test_get_edges_in_groups.jl")
    
    # Only run performance tests if BenchmarkTools is available
    if RUN_BENCHMARKS
        @testset "Performance Regression" begin
            include("test_performance_regression.jl")
        end
    end
end
