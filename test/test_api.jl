@testset "test api" begin
    using Distributions
    A = rand(0:1, 40, 40)
    for i in 1:40
        A[i, i] = 0
    end

    g = Observations(Symmetric(A), Uniform(-1, 1))
    sbm_fitted = nethist(g; h = 10, max_iter = 10)

    @test eltype(sbm_fitted) == typeof(Uniform(-1, 1))
    @test size(sbm_fitted) == (4,4)

    sbm_discretised, discretizer = nethist_discretised(
        g; number_levels = 5, h = 10, max_iter = 10)
    @test eltype(sbm_discretised) == typeof(Categorical(5))
    @test size(sbm_discretised) == (4,4)
end
