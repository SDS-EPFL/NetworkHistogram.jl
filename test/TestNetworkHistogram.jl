module TestNetworkHistogram

import NetworkHistogram as NH
using Test

function to_default_assignment(a_specialised::NH.Assignment{T, B}) where {T, B}
    return NH.Assignment(a_specialised.group_size, a_specialised.node_labels)
end

to_default_assignment(a::NH.Assignment{T, Nothing}) where {T} = a

function test_swap_revertible(
        a::NH.Assignment, swap::NH.Swap, g::NH.Observations)
    a_test = deepcopy(a)
    NH.apply_swap!(a_test, swap)
    @test NH.get_group_of_vertex(a, swap.index1) ==
          NH.get_group_of_vertex(a_test, swap.index2)
    @test NH.get_group_of_vertex(a, swap.index2) ==
          NH.get_group_of_vertex(a_test, swap.index1)

    # force recomputation of the log likelihood using default assignment
    a_new = to_default_assignment(a_test)
    @test NH.log_likelihood(a_new, g) ≈ NH.log_likelihood(a_test, g)

    # revert the swap and check if the assignment is the same as before
    NH.revert_swap!(a_test, swap)
    @test a == a_test
    @test NH.log_likelihood(a, g) ≈ NH.log_likelihood(a_test, g)
end

end
