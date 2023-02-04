@testset "Pipeline" begin
    @testset "dummy run" begin
        A = [0 0 1 0 1 0 1 1 0 1
             0 0 1 1 1 1 1 1 0 0
             1 1 0 1 0 0 0 0 1 0
             0 1 1 0 1 0 1 0 0 0
             1 1 0 1 0 0 1 0 0 1
             0 1 0 0 0 0 0 1 0 0
             1 1 0 1 1 0 0 1 0 1
             1 1 0 0 0 1 1 0 0 1
             0 0 1 0 0 0 0 0 0 1
             1 0 0 0 1 0 1 1 1 0]
        graphist = graphhist(A; h = 0.5)
    end

    @testset "associative stochastic block model" begin
        adjacencies = load(pwd() * "/test_files/sbm.jld")

        for (name, adjacency) in adjacencies
            @testset "$name" begin
                graphist, history = graphhist(adjacency; h = 0.3,
                                              stop_rule = PreviousBestValue(100),
                                              starting_assignment_rule = OrderedStart())
                @test all(graphist.θ .>= 0.0)
            end
        end
    end
end
