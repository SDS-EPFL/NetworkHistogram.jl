abstract type AcceptRule end
struct Strict <: AcceptRule end

"""
    accept_reject_update!(a::Assignment, swap::Swap, g, accept_rule::AcceptRule)


Perform the swap and accept it if it improves the likelihood of the assignment. `a` will
be updated in place if the swap is accepted.

# Implemented rules
- `Strict()`: Accept the proposal if it has a higher likelihood than the current assignment.
"""
accept_reject_update!

function accept_reject_update!(a::Assignment, swap::Swap, ::Strict)
    current_score = loglikelihood(a)
    apply_swap!(a, swap)
    if loglikelihood(a) <= current_score
        revert_swap!(a, swap)
    end
    return nothing
end
