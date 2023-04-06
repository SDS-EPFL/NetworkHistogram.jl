@testset "Error handling" begin
    @testset "Adjacency matrix" begin
        As = [
            [0 1
             0 0], [1 1
                    1 0], [0 2
                           2 0], [0 1
                                  1 0
                                  0 1],
        ]
        for A in As
            @test_throws AssertionError graphhist(A, h = 2)
            @test_throws AssertionError graphhist(Bool.(min.(A, 1)), h = 2)
            @test_throws AssertionError graphhist(Float64.(A), h = 2)
        end
        @test_throws AssertionError graphhist(["0" "1"; "1" "0"], h = 2)
    end
    @testset "maxitr" begin
        @test_throws AssertionError graphhist([0 1; 1 0], h = 2,
                                              maxitr = -1)
    end
    @testset "h" begin for h in (3, -1, 1.1, -0.1)
        @test_throws AssertionError graphhist([0 1; 1 0], h = h)
    end end
end
