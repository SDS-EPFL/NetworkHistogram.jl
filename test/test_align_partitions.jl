@testset "align_partitions" begin
    @testset "Basic alignment" begin
        x = [1, 1, 2, 2, 3]
        y = [2, 2, 1, 1, 3]

        y_aligned, mapping = align_partitions(x, y)

        @test length(y_aligned) == length(y)
        @test y_aligned == [1, 1, 2, 2, 3]
        @test mapping[2] == 1
        @test mapping[1] == 2
        @test mapping[3] == 3
    end

    @testset "Perfect match" begin
        x = [1, 1, 2, 2, 3, 3]
        y = [1, 1, 2, 2, 3, 3]

        y_aligned, mapping = align_partitions(x, y)

        @test y_aligned == x
        @test mapping == Dict(1 => 1, 2 => 2, 3 => 3)
    end

    @testset "Complete permutation" begin
        x = [1, 1, 2, 2, 3, 3]
        y = [3, 3, 1, 1, 2, 2]

        y_aligned, mapping = align_partitions(x, y)

        @test y_aligned == x
        @test mapping == Dict(3 => 1, 1 => 2, 2 => 3)
    end

    @testset "Different numbers of clusters" begin
        # x has 3 clusters, y has 2
        x = [1, 1, 2, 2, 3, 3]
        y = [1, 1, 1, 1, 2, 2]

        y_aligned, mapping = align_partitions(x, y)

        @test length(y_aligned) == length(y)
        # The alignment should maximize overlap
        # y cluster 1 should map to x cluster with most overlap (1 or 2)
        # y cluster 2 should map to x cluster 3
        @test sum(y_aligned[1:4] .== x[1:4]) ≥ 2  # Good overlap for first 4
    end

    @testset "Single cluster" begin
        x = [1, 1, 1, 1]
        y = [1, 1, 1, 1]

        y_aligned, mapping = align_partitions(x, y)

        @test y_aligned == x
        @test mapping == Dict(1 => 1)
    end

    @testset "Non-contiguous labels" begin
        x = [10, 10, 20, 20, 30, 30]
        y = [5, 5, 15, 15, 25, 25]

        y_aligned, mapping = align_partitions(x, y)

        @test length(y_aligned) == length(y)
        # Should create optimal matching
        @test all(l -> l ∈ [10, 20, 30], y_aligned)
    end

    @testset "Error handling" begin
        x = [1, 2, 3]
        y = [1, 2]

        @test_throws Exception align_partitions(x, y)
    end

    @testset "Preserves partition structure" begin
        x = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        y = [2, 2, 2, 3, 3, 3, 1, 1, 1]

        y_aligned, mapping = align_partitions(x, y)

        # Check that elements in same cluster in y stay together in y_aligned
        @test length(unique(y_aligned[1:3])) == 1
        @test length(unique(y_aligned[4:6])) == 1
        @test length(unique(y_aligned[7:9])) == 1

        # Should achieve perfect alignment after relabeling
        @test y_aligned == x
    end
end
