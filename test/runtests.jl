using Test
using LinearAlgebra, SparseArrays
using NetworkHistogram
@testset "Tests" begin
    @testset "test can run" begin
        @test 1 == 1
    end


    @testset "Edge list tests" begin
        A = Symmetric(sprand(20,20,0.5))
        edgelist = EdgeList(A)

        for j in 1:20
            for i in 1:20
                if A[i,j] != 0
                    nv_j, val_j = neighbors(edgelist, j)
                    @test i in nv_j
                    @test A[i,j] == val_j[findfirst(x -> x == i, nv_j)]
                end
            end
        end

        @test eltype(edgelist) == eltype(A)
        @test nodes(edgelist) == size(A,1)
    end

end
