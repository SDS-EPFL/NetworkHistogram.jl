@testset begin
    "Pipeline"
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
    graphist = NetworkHistogram.graphhist(A; h = 0.5)
end

@testset begin
    "SBM"
    adjacencies = load(pwd() * "/test_files/sbm.jld")

    for (name, adjacency) in adjacencies
        @testset begin
            name
            graphist,_ = NetworkHistogram.graphhist(adjacency; h = 0.3,
                                                  stop_rule = NetworkHistogram.PreviousBestValue(10))
            println(graphist.θ)
            @test all(graphist.θ .>= 0.0)
        end
    end
end
