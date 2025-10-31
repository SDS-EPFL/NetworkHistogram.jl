using Test
using LinearAlgebra, SparseArrays
using NetworkHistogram

@testset "Tests" begin
    include("test_symarray.jl")
    include("test_pseudo_suff_stats.jl")
    include("test_hist_dist.jl")
    include("test_align_partitions.jl")
end
