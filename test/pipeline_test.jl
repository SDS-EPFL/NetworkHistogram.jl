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
        @testset "run bandwidth float" begin
            estimated = graphhist(A; h = 0.5)
            @test all(estimated.graphhist.θ .>= 0.0)
            @test all(estimated.graphhist.θ .<= 1.0)
            @test size(estimated.graphhist.θ) == (2, 2)
        end
        @testset "run bandwidth int" begin
            estimated = graphhist(A; h = 5)
            @test all(estimated.graphhist.θ .>= 0.0)
            @test all(estimated.graphhist.θ .<= 1.0)
            @test size(estimated.graphhist.θ) == (2, 2)
        end
    end

    @testset "associative stochastic block model" begin
        adjacencies = load(pwd() * "/test_files/sbm.jld")

        for (name, adjacency) in adjacencies
            @testset "$name" begin
                estimated, history = graphhist(adjacency; h = 0.3,
                                               stop_rule = PreviousBestValue(100),
                                               starting_assignment_rule = OrderedStart())
                @test all(estimated.θ .>= 0.0)
                estimated, history = graphhist(adjacency; h = 0.3,
                                               stop_rule = PreviousBestValue(100),
                                               starting_assignment_rule = OrderedStart(),
                                               record_trace = false)
                @test all(estimated.θ .>= 0.0)
            end
        end
    end
end
