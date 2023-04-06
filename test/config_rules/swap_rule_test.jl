import NetworkHistogram: select_swap
@testset "swap rule" begin
    A, node_labels, group_size, assignment = make_simple_example()
    x = select_swap(assignment, A, RandomNodeSwap())
    @test x isa Tuple{Int,Int}
    @test all(1 .≤ x .≤ size(A,1))
end
