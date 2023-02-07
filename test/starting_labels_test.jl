

"""
    test_basic_node_labels(node_labels, group_size)

Test that the node labels are valid:
    - correct number of labels
    - labels are positive
    - labels are within the range of number of groups
"""
function test_basic_node_labels(node_labels, group_size)
    @test length(node_labels) == sum(group_size)
    @test all(node_labels .> 0)
    @test all(node_labels .<= length(group_size))
    for (i, group_s) in enumerate(group_size)
        @test count(x -> x == i, node_labels) == group_s
    end
end

@testset "Initial node labels" begin

    A = [0 1 1 1 0 0 1 0
         1 0 1 1 0 0 0 0
         1 1 0 0 0 0 0 0
         1 1 0 0 0 0 0 1
         0 0 0 0 0 1 1 1
         0 0 0 0 1 0 1 1
         1 0 0 0 1 1 0 0
         0 0 0 1 1 1 0 0]
    h = 0.5

    @testset "random start" begin
        node_labels, group_size = NetworkHistogram.initialise_node_labels(A, h, NetworkHistogram.RandomStart())
        test_basic_node_labels(node_labels, group_size)
    end

    @testset "ordered start" begin
        node_labels, group_size = NetworkHistogram.initialise_node_labels(A, h, NetworkHistogram.OrderedStart())
        test_basic_node_labels(node_labels, group_size)
        @test node_labels == sort(node_labels)
    end

    @testset "eigenvalue start" begin
        node_labels, group_size = NetworkHistogram.initialise_node_labels(A, h, NetworkHistogram.EigenStart())
        test_basic_node_labels(node_labels, group_size)
    end

end