abstract type AcceptRule end
struct Strict <: AcceptRule end

"""
    accept_reject_update!(history::GraphOptimizationHistory, iteration::Int,
                          proposal::Assignment,
                          current::Assignment, accept_rule::AcceptRule)


Return the updated `current` assignment based on the `accept_rule`.

# Implemented rules
- `Strict()`: Accept the proposal if it has a higher likelihood than the current assignment.
"""
accept_reject_update!

function accept_reject_update!(history::GraphOptimizationHistory, iteration::Int,
        proposal::Assignment,
        current::Assignment, ::Strict)
    if proposal.likelihood > current.likelihood
        deepcopy!(current, proposal)
    end

    update_current!(history, iteration, current.likelihood)
    return current
end


function accept_reject_update!(history::GraphOptimizationHistory, iteration::Int,
        proposal::Assignments,
        current::Assignments, ::Strict)
    if proposal.likelihood > current.likelihood
        deepcopy!(current, proposal)
    end

    update_current!(history, iteration, current.likelihood)
    return current
end
