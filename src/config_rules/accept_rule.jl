abstract type AcceptRule end
struct Strict <: AcceptRule end

function accept_reject_update!(history::MVHistory, iteration::Int, proposal::Assignment,
                               current::Assignment, accept_rule::Strict)
    if proposal.likelihood[1] > current.likelihood[1]
        deepcopy!(current, proposal)
    end

    push!(history, :current_likelihood, iteration, current.likelihood[1])
    return current
end
