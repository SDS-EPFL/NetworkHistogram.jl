import NetworkHistogram: stopping_rule
@testset "stop rule" begin
    A, node_labels, group_size, proposal = make_simple_example()

    proposal.likelihood = -0.1
    current = deepcopy(proposal)
    best = deepcopy(proposal)

    histories = [
        initialize_history(best, current, proposal, Val{true}()),
        initialize_history(best, current, proposal, Val{false}()),
    ]
    for history in histories
        for i in 1:4
            NetworkHistogram.update_current!(history, i, 0.0)
        end
        @test stopping_rule(history, PreviousBestValue(3)) == true
        @test stopping_rule(history, PreviousBestValue(4)) == false
    end
end
