abstract type AcceptRule end
struct Strict <: AcceptRule end

function accept_reject_update!(history::GraphOptimizationHistory, iteration::Int,
                               proposal::Assignment,
                               current::Assignment, accept_rule::Strict)
    if proposal.likelihood > current.likelihood
        deepcopy!(current, proposal)
    end

    update_current!(history, iteration, current.likelihood)
    return current
end
