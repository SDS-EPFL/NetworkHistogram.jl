import NetworkHistogram: accept_reject_update!
@testset "accept rule" begin
    accept_rules = [Strict()]

    iteration = 3
    proposal = NetworkHistogram.Assignment([0 1; 1 0], [1, 2], 1)
    proposal.likelihood = 0.0
    current = deepcopy(proposal)
    best = deepcopy(proposal)

    histories = [
        initialise_history(best, current, proposal, Val{true}()),
        initialise_history(best, current, proposal, Val{false}()),
    ]
    proposal_likelihoods = [-0.1, 0.1]
    for history in histories, accept_rule in accept_rules, l in proposal.likelihoods
        proposal.likelihood = l
        accept_reject_update!(history, iteration, proposal, current, accept_rule)
        @test current.likelihood == max(l, 0.0)
        current.likelihoood = 0.0
    end
end
