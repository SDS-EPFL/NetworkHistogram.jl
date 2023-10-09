@testset "multilayer" begin
    @testset "test initialisation of assignments" begin
        A, labels, group_size, assignment = make_simple_example()
        A2, _,_, assignment2 = make_multivariate_example()
        @test all(assignment.estimated_theta .== assignment2.estimated_theta[:,:,2])
        @test all(assignment.realized .== assignment2.realized[:, :, 2])
        @test assignment.likelihood == assignment2.likelihood
        @test sum(assignment2.estimated_theta) ≈ size(assignment2.estimated_theta,1)^2
    end
end
