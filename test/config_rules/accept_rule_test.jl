import NetworkHistogram: accept_reject_update!, initialize_history
@testset "accept rule" begin

    iteration = 3
    A, node_labels, group_size, proposal = make_simple_example()

    proposal.likelihood = 0.0
    current = deepcopy(proposal)
    best = deepcopy(proposal)

    histories = [
        initialize_history(best, current, proposal, Val{true}()),
        initialize_history(best, current, proposal, Val{false}()),
    ]
    test_likelihoods = [-0.1, 0.1]
    for history in histories, lik in test_likelihoods
        proposal.likelihood = lik # set proposal
        accept_reject_update!(history, iteration, proposal, current, Strict())

        @testset "Strict with history is $(typeof(history).name.name), likelihood=$lik" begin
            @test current.likelihood == max(lik, 0.0) # should have accepted if better
            if history isa NetworkHistogram.TraceHistory
                @test get(history.history, :current_likelihood)[1][end] == iteration
                @test get(history.history, :current_likelihood)[2][end] == current.likelihood
            else
                @test history.current_iteration == iteration
            end
        end

        current.likelihood = 0.0 # reset for next example
        iteration += 1 # otherwise will get an error from history
    end
end
