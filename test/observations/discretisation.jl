using NetworkHistogram

@testset "discretisation" begin
    using Distributions
    A = rand(-1:1, 20, 20)
    for i in 1:20
        A[i, i] = 0
    end
    g, discretizer = Observations(A,Uniform(-1,1))
    discretised_g = discretise(g; number_levels = 5)
    @test size(discretised_g.graph) == size(g.graph)
    @test discretised_g.dist_ref == Categorical(6)
    @test all(discretised_g.graph .∈ Ref(0:5))
end
