@testset "swap test" begin
    # given a::Assignment and s::Swap
    # save copy assignment
    a_ref = deepcopy(a)
    revert_swap!(make_swap!(a), s)
    @assert a == a_ref
end
