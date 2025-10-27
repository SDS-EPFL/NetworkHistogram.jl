using Test
using NetworkHistogram
using SparseArrays
using LinearAlgebra
using StaticArrays

@testset "SymArray Array Interface" begin
    @testset "Construction and basic properties" begin
        # Test construction with scalar
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 1.0)
        @test a isa AbstractArray{Float64, 2}
        @test size(a) == (3, 3)
        @test length(a) == 9
        @test axes(a) == (1:3, 1:3)
        @test eltype(a) == Float64

        # Test construction with zeros
        b = SymArray{Float64}(undef, 5, 5)
        fill!(b, 0.0)
        @test size(b) == (5, 5)
        @test all(b[i, j] == 0.0 for i in 1:5 for j in 1:5)

        # Test dimension validation
        @test_throws ArgumentError SymArray{Float64}(undef, 3, 4)
    end

    @testset "Indexing - getindex and setindex!" begin
        a = SymArray{Float64}(undef, 4, 4)
        fill!(a, 0.0)

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
        a = SymArray{Float64}(undef, 5, 5)
        fill!(a, 0.0)

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
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 5.0)

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
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 0.0)
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
        d = SymArray{Float64}(undef, 4, 4)
        fill!(d, 0.0)
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
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 2.0)

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
        b = SymArray{Float64}(undef, 3, 3)
        fill!(b, 0.0)
        b[1, 1] = 5.0
        b[2, 3] = -3.0
        @test maximum(b) == 5.0
        @test minimum(b) == -3.0
    end

    @testset "Mathematical operations" begin
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 2.0)
        b = SymArray{Float64}(undef, 3, 3)
        fill!(b, 3.0)

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
        g = SymArray{Float64}(undef, 3, 3)
        fill!(g, -2.0)
        h = abs.(g)
        @test h isa SymArray
        @test all(h[i, j] == 2.0 for i in 1:3, j in 1:3)

        # Test with mixed values
        m = SymArray{Float64}(undef, 3, 3)
        fill!(m, 0.0)
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
        p = SymArray{Float64}(undef, 3, 3)
        fill!(p, 0.0)
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
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 1.0)
        # Only upper triangle is stored: 6 elements
        # [1,1], [1,2], [1,3], [2,2], [2,3], [3,3]
        @test sum_tri_with_diag(a) == 6.0

        b = SymArray{Float64}(undef, 4, 4)
        fill!(b, 2.0)
        # Upper triangle has 10 elements for 4x4
        @test sum_tri_with_diag(b) == 20.0

        # Verify it's different from full sum (which counts off-diag twice)
        # Full sum would be 2*n*(n-1)/2 + n for value v
        # = v*(n^2-n+n) = v*n^2
        # While sum_tri_with_diag gives v*n*(n+1)/2
    end

    @testset "Type stability" begin
        # Float64
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 1.0)
        @test typeof(a[1, 1]) == Float64

        # Int
        b = SymArray{Int}(undef, 3, 3)
        fill!(b, 1)
        @test typeof(b[1, 1]) == Int

        # Float32
        c = SymArray{Float32}(undef, 3, 3)
        fill!(c, 1.0f0)
        @test typeof(c[1, 1]) == Float32
    end

    @testset "Sparse matrix properties" begin
        a = SymArray{Float64}(undef, 10, 10)
        fill!(a, 0.0)
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
        a = SymArray{Float64}(undef, 1, 1)
        fill!(a, 5.0)
        @test size(a) == (1, 1)
        @test a[1, 1] == 5.0
        a[1, 1] = 10.0
        @test a[1, 1] == 10.0

        # Large diagonal
        b = SymArray{Float64}(undef, 100, 100)
        fill!(b, 0.0)
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
        sin_a_bis = sin.(a)
        for sin_test in (sin_a, sin_a_bis)
            @test sin_test isa SymArray
            @test all(sin_test[i, j] == sin(2.0) for i in 1:3, j in 1:3)
        end
    end

    @testset "Broadcasting with regular arrays" begin
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 2.0)
        M = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        # SymArray + Matrix should return Matrix (follows Matrix type)
        result1 = a .+ M
        @test result1 isa Matrix{Float64}
        @test !(result1 isa SymArray)

        # Matrix + SymArray should also return Matrix
        result2 = M .+ a
        @test result2 isa Matrix{Float64}
        @test !(result2 isa SymArray)

        # Check values are correct
        for i in 1:3, j in 1:3
            @test result1[i, j] ≈ 2.0 + M[i, j]
            @test result2[i, j] ≈ M[i, j] + 2.0
        end

        # SymArray + scalar should still return SymArray
        result3 = a .+ 5.0
        @test result3 isa SymArray

        # SymArray + SymArray should return SymArray
        b = make_sym_init(3, 3.0)
        result4 = a .+ b
        @test result4 isa SymArray
    end

    @testset "SymArray broadcast with Matrix returns Matrix" begin
        # Create a SymArray and a regular Matrix
        a = SymArray{Float64}(undef, 3, 3)
        fill!(a, 2.0)
        M = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        # SymArray + Matrix should return Matrix
        result1 = a .+ M
        @test result1 isa Matrix{Float64}
        @test !(result1 isa SymArray)
        @test size(result1) == (3, 3)

        # Matrix + SymArray should also return Matrix
        result2 = M .+ a
        @test result2 isa Matrix{Float64}
        @test !(result2 isa SymArray)

        # Check values are correct
        for i in 1:3, j in 1:3
            @test result1[i, j] ≈ 2.0 + M[i, j]
            @test result2[i, j] ≈ M[i, j] + 2.0
        end

        # SymArray + scalar should still return SymArray
        result3 = a .+ 5.0
        @test result3 isa SymArray
        @test all(result3[i, j] ≈ 7.0 for i in 1:3, j in 1:3)

        # SymArray + SymArray should return SymArray
        b = SymArray{Float64}(undef, 3, 3)
        fill!(b, 3.0)
        result4 = a .+ b
        @test result4 isa SymArray
        @test all(result4[i, j] ≈ 5.0 for i in 1:3, j in 1:3)

        # Chained operations with scalars should still work
        result5 = (a .+ 1) .* 2
        @test result5 isa SymArray
        @test all(result5[i, j] ≈ 6.0 for i in 1:3, j in 1:3)

        a_ones = SymArray{Float64}(undef, 3, 3)
        fill!(a_ones, 1.0)
        result_sum_two_matrices = a_ones .+ M .+ M
        @test result_sum_two_matrices isa Matrix{Float64}
        @test all(result_sum_two_matrices[i, j] ≈ 1 + 2 * M[i, j] for i in 1:3, j in 1:3)
    end
end
