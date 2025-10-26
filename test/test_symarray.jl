using Test
using NetworkHistogram
using SparseArrays
using LinearAlgebra
using StaticArrays

@testset "SymArray Array Interface" begin
    @testset "Construction and basic properties" begin
        # Test construction with scalar
        a = SymArray(3, 1.0)
        @test a isa AbstractArray{Float64, 2}
        @test size(a) == (3, 3)
        @test length(a) == 9
        @test axes(a) == (1:3, 1:3)
        @test eltype(a) == Float64

        # Test construction with zeros
        b = SymArray(5, 0.0)
        @test size(b) == (5, 5)
        @test all(b[i, j] == 0.0 for i in 1:5 for j in 1:5)

        # Test dimension validation
        @test_throws ArgumentError SymArray(0, 1.0)
        @test_throws ArgumentError SymArray(-1, 1.0)
    end

    @testset "Indexing - getindex and setindex!" begin
        a = SymArray(4, 0.0)

        # Test setindex! in upper triangle
        a[1, 2] = 5.0
        @test a[1, 2] == 5.0
        @test a[2, 1] == 5.0  # Symmetry

        # Test setindex! in lower triangle (should set upper)
        a[3, 2] = 7.0
        @test a[2, 3] == 7.0
        @test a[3, 2] == 7.0

        # Test diagonal
        a[2, 2] = 3.0
        @test a[2, 2] == 3.0

        # Test bounds checking
        @test_throws BoundsError a[0, 1]
        @test_throws BoundsError a[5, 1]
        @test_throws BoundsError a[1, 5]
    end

    @testset "Symmetry property" begin
        a = SymArray(5, 0.0)

        # Set values and verify symmetry
        for i in 1:5
            for j in 1:5
                val = i * 10 + j
                a[i, j] = val
                @test a[i, j] == a[j, i]
            end
        end
    end

    @testset "Construction from matrix" begin
        # Test from symmetric matrix
        M = [1.0 2.0 3.0;
             2.0 4.0 5.0;
             3.0 5.0 6.0]
        a = SymArray(M)

        @test size(a) == (3, 3)
        for i in 1:3, j in 1:3
            @test a[i, j] == M[i, j]
        end

        # Test non-square matrix throws error
        @test_throws ArgumentError SymArray([1.0 2.0; 3.0 4.0; 5.0 6.0])
    end

    @testset "convert functions" begin
        # Test conversion to SymArray
        M = [1.0 2.0; 2.0 4.0]
        a = convert(SymArray{Float64}, M)
        @test a isa SymArray{Float64}
        @test a[1, 1] == 1.0
        @test a[1, 2] == 2.0
        @test a[2, 2] == 4.0

        # Test conversion to AbstractMatrix
        b = convert(Matrix{Float64}, a)
        @test b isa Matrix{Float64}
        @test b == M
        @test b[1, 2] == b[2, 1]  # Verify symmetry
    end

    @testset "similar function" begin
        a = SymArray(3, 5.0)

        # Test similar without type
        b = similar(a)
        @test size(b) == size(a)
        @test eltype(b) == eltype(a)
        @test b isa SymArray{Float64}

        # Test similar with type
        c = similar(a, Int)
        @test size(c) == size(a)
        @test eltype(c) == Int
        @test c isa SymArray{Int}

        # Test similar with type and dimensions
        d = similar(a, Float32, (4, 4))
        @test size(d) == (4, 4)
        @test eltype(d) == Float32

        # Test non-square dimensions throw error
        @test_throws ArgumentError similar(a, Float64, (3, 4))
    end

    @testset "copy! and deepcopy!" begin
        a = SymArray(3, 0.0)
        a[1, 1] = 1.0
        a[1, 2] = 2.0
        a[2, 3] = 5.0

        b = similar(a)
        copy!(b, a)

        @test b[1, 1] == 1.0
        @test b[1, 2] == 2.0
        @test b[2, 1] == 2.0
        @test b[2, 3] == 5.0
        @test b[3, 2] == 5.0

        # Test dimension mismatch
        d = SymArray(4, 0.0)
        @test_throws DimensionMismatch copy!(d, a)

        # Test deepcopy!
        src = SymArray{Vector{Int}}(undef, 4, 4)
        for j in 1:4, i in j:4
            src[i, j] = [i, j]
        end

        # on unassigned dest
        dest = similar(src)
        deepcopy!(dest, src)
        for j in 1:4, i in j:4
            @test dest[i, j] == src[i, j]
            @test !(dest[i, j] === src[i, j])  # Ensure deep copy
        end

        # on assigned dest
        dest2 = similar(src)
        for j in 1:4, i in j:4
            dest2[i, j] = [-1, -1]
        end
        deepcopy!(dest2, src)
        for j in 1:4, i in j:4
            @test dest2[i, j] == src[i, j]
            @test !(dest2[i, j] === src[i, j])  # Ensure deep copy
        end
    end

    @testset "Array operations" begin
        a = SymArray(3, 2.0)

        # Test iteration
        count = 0
        for val in a
            @test val == 2.0
            count += 1
        end
        @test count == 9

        # Test sum
        @test sum(a) == 18.0

        # Test all/any
        @test all(x -> x == 2.0, a)
        @test any(x -> x == 2.0, a)

        # Test maximum/minimum
        b = SymArray(3, 0.0)
        b[1, 1] = 5.0
        b[2, 3] = -3.0
        @test maximum(b) == 5.0
        @test minimum(b) == -3.0
    end

    @testset "Mathematical operations" begin
        a = SymArray(3, 2.0)
        b = SymArray(3, 3.0)

        # Element-wise operations (using broadcasting)
        c = a .+ b
        @test c isa SymArray
        @test all(c[i, j] == 5.0 for i in 1:3, j in 1:3)

        d = a .* 2
        @test d isa SymArray
        @test all(d[i, j] == 4.0 for i in 1:3, j in 1:3)

        # Test subtraction
        e = b .- a
        @test e isa SymArray
        @test all(e[i, j] == 1.0 for i in 1:3, j in 1:3)

        # Test division
        f = b ./ 2.0
        @test f isa SymArray
        @test all(f[i, j] == 1.5 for i in 1:3, j in 1:3)

        # Test unary operations
        g = SymArray(3, -2.0)
        h = abs.(g)
        @test h isa SymArray
        @test all(h[i, j] == 2.0 for i in 1:3, j in 1:3)

        # Test with mixed values
        m = SymArray(3, 0.0)
        m[1, 1] = 1.0
        m[1, 2] = 2.0
        m[2, 2] = 3.0
        m[1, 3] = 4.0
        m[2, 3] = 5.0
        m[3, 3] = 6.0

        n = m .+ 10.0
        @test n isa SymArray
        @test n[1, 1] == 11.0
        @test n[1, 2] == 12.0
        @test n[2, 1] == 12.0  # Symmetry
        @test n[2, 2] == 13.0
        @test n[3, 3] == 16.0

        # Test operations between two SymArrays with different values
        p = SymArray(3, 0.0)
        p[1, 1] = 10.0
        p[2, 2] = 20.0
        p[3, 3] = 30.0

        q = m .+ p
        @test q isa SymArray
        @test q[1, 1] == 11.0
        @test q[2, 2] == 23.0
        @test q[3, 3] == 36.0
        @test q[1, 2] == 2.0
        @test q[2, 1] == 2.0
    end

    @testset "Special case: sum_tri_with_diag" begin
        a = SymArray(3, 1.0)
        # Only upper triangle is stored: 6 elements
        # [1,1], [1,2], [1,3], [2,2], [2,3], [3,3]
        @test sum_tri_with_diag(a) == 6.0

        b = SymArray(4, 2.0)
        # Upper triangle has 10 elements for 4x4
        @test sum_tri_with_diag(b) == 20.0

        # Verify it's different from full sum (which counts off-diag twice)
        # Full sum would be 2*n*(n-1)/2 + n for value v
        # = v*(n^2-n+n) = v*n^2
        # While sum_tri_with_diag gives v*n*(n+1)/2
    end

    @testset "Type stability" begin
        # Float64
        a = SymArray(3, 1.0)
        @test typeof(a[1, 1]) == Float64

        # Int
        b = SymArray(3, 1)
        @test typeof(b[1, 1]) == Int

        # Float32
        c = SymArray(3, 1.0f0)
        @test typeof(c[1, 1]) == Float32
    end

    @testset "Sparse matrix properties" begin
        a = SymArray(10, 0.0)
        # Initially all elements are stored (including zeros)
        # Set only a few elements to non-zero
        a[1, 5] = 3.0
        a[3, 7] = 4.0
        a[9, 9] = 5.0

        # Verify values are correct (symmetry)
        @test a[1, 5] == 3.0
        @test a[5, 1] == 3.0
        @test a[3, 7] == 4.0
        @test a[7, 3] == 4.0
        @test a[9, 9] == 5.0
        @test a[2, 2] == 0.0
    end

    @testset "Edge cases" begin
        # 1x1 matrix
        a = SymArray(1, 5.0)
        @test size(a) == (1, 1)
        @test a[1, 1] == 5.0
        a[1, 1] = 10.0
        @test a[1, 1] == 10.0

        # Large diagonal
        b = SymArray(100, 0.0)
        for i in 1:100
            b[i, i] = Float64(i)
        end
        @test b[50, 50] == 50.0
        @test b[99, 99] == 99.0
    end

    @testset "Broadcasting" begin
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 2.0)
        b = @. a + 2.0
        @test b isa SymArray
        @test all(b[i, j] == 4.0 for i in 1:3, j in 1:3)

        c = b ./ a
        @test c isa SymArray
        @test all(c[i, j] == 2.0 for i in 1:3, j in 1:3)

        sin_a = @. sin(a)
        @test sin_a isa SymArray
        @test all(sin_a[i, j] == sin(2.0) for i in 1:3, j in 1:3)
    end
end
