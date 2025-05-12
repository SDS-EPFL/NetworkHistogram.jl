using Test
using NetworkHistogram

@testset "get_edges_in_groups behavior" begin
    # Simple 4-node undirected graph
    # 1-2, 1-3, 2-4, 3-4
    A = [0 1 1 0;
         1 0 0 1;
         1 0 0 1;
         0 1 1 0]
    edgelist = NetworkHistogram.EdgeList(A)
    node_labels = [1, 1, 2, 2]  # nodes 1,2 in group 1; 3,4 in group 2

    # Test within-group edges (group 1)
    edges_1_1 = NetworkHistogram.get_edges_in_groups(node_labels, edgelist, 1, 1)
    @test length(edges_1_1) == 1  # Only edge (1,2)
    @test edges_1_1[1] == 1  # A[1,2] == 1

    # Test within-group edges (group 2)
    edges_2_2 = NetworkHistogram.get_edges_in_groups(node_labels, edgelist, 2, 2)
    @test length(edges_2_2) == 1  # Only edge (3,4)
    @test edges_2_2[1] == 1  # A[3,4] == 1

    # Test between-group edges (1,2)
    edges_1_2 = NetworkHistogram.get_edges_in_groups(node_labels, edgelist, 1, 2)
    # Edges: (1,3), (2,4)
    @test length(edges_1_2) == 4
    @test sort(edges_1_2) == [0, 0, 1, 1]  # Both edges exist

    # Test symmetry: get_edges_in_groups(2,1) == get_edges_in_groups(1,2)
    edges_2_1 = NetworkHistogram.get_edges_in_groups(node_labels, edgelist, 2, 1)
    @test sort(edges_2_1) == sort(edges_1_2)
end
