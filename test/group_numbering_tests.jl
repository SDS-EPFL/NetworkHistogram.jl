@testset "Group Numbers" begin
    A, _, _, _ = make_simple_example()
    h = 0.5
    node_labels, group_size = NetworkHistogram.initialize_node_labels(A, h,
        RandomStart())
    group_size_new = NetworkHistogram.GroupSize(node_labels)
    @test size(group_size_new) == size(group_size)
    @test group_size_new.group_number == group_size.group_number
    @test group_size_new.number_groups == group_size.number_groups
end
