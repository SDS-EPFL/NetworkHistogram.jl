@testset "Data utils" begin
    A = [0 0 0 1
         0 0 0 0
         1 0 0 1
         0 0 1 0]

    @testset "drop isolated vertices" begin
        B = NetworkHistogram.drop_isolated_vertices(A)
        @test B == [0 1 1
                    1 0 1
                    0 1 0]
    end
end
