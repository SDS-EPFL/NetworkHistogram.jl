@testset "Pipeline" begin
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
    @testset "dummy run" begin
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
        @testset "run with automatic bandwidth" begin
            estimated = graphhist(A)
            @test all(estimated.graphhist.θ .>= 0.0)
            @test all(estimated.graphhist.θ .<= 1.0)
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

    @testset "multilayer run" begin
        @testset "2 layers perfectly correlated" begin
            A_2 = cat(A, A, dims = 3)
            estimated, history = graphhist(A_2; h = 0.5)
            @test all(estimated.θ .>= 0.0)
            @test all(estimated.θ .<= 1.0)
            @test size(estimated.θ) == (2, 2, 4)
        end
        @testset "run with automatic bandwidth" begin
            A_2 = cat(A, A, dims = 3)
            estimated, history = graphhist(A_2)
            @test all(estimated.θ .>= 0.0)
            @test all(estimated.θ .<= 1.0)
        end

        @testset "2 layers perfectly anti-correlated" begin
            A_2 = cat(A, abs.(A .- 1), dims = 3)
            for i in 1:size(A, 1)
                A_2[i, i, 2] = 0
            end
            estimated, history = graphhist(A_2; h = 0.5)
            @test all(estimated.θ .>= 0.0)
            @test all(estimated.θ .<= 1.0)
            @test size(estimated.θ) == (2, 2, 4)
        end

        @testset "3 layers" begin
            A_3 = cat(A, A, A, dims = 3)
            estimated, history = graphhist(A_3; h = 0.5)
            @test all(estimated.θ .>= 0.0)
            @test all(estimated.θ .<= 1.0)
            @test size(estimated.θ) == (2, 2, 8)
        end
    end
end
