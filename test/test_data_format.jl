@testset "Edge list tests" begin
    using Random
    Random.seed!(1234)
    A = Symmetric(sprand(20,20,0.5))
    edgelist = EdgeList(A)

    for j in 1:20
        nv_j, val_j = neighbors(edgelist, j)
        for i in 1:20
            if A[i,j] == 0
                @test i ∉ nv_j
            end
            if A[i,j] != 0
                @test i in nv_j
                @test A[i,j] == val_j[findfirst(x -> x == i, nv_j)]
            end
        end
    end

    @test eltype(edgelist) == eltype(A)
    @test nodes(edgelist) == size(A,1)
end
