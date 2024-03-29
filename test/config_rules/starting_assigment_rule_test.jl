import NetworkHistogram: initialize_node_labels

@testset "starting assignment rules" begin
    @testset "starting assigment rule simple graphs" begin
        A, _, _, _ = make_simple_example()
        for method in (OrderedStart(),
            RandomStart(),
            EigenStart(),
            DistStart())
            node_labels, group_size = initialize_node_labels(A, 4, method)
            if method isa OrderedStart
                @test sort(node_labels) == node_labels
            end
            @test all(sum(n -> n == j, node_labels) == group_size[j]
                      for j in unique(node_labels))
        end
    end

    @testset "starting assigment rule multilayer graphs" begin
        A, _, _, _ = make_multivariate_example()
        for method in (OrderedStart(),
            RandomStart(),
            EigenStart(),
            DistStart())
            node_labels, group_size = initialize_node_labels(A, 4, method)
            if method isa OrderedStart
                @test sort(node_labels) == node_labels
            end
            @test all(sum(n -> n == j, node_labels) == group_size[j]
                      for j in unique(node_labels))
        end
    end
end
