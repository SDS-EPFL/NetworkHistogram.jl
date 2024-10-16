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

function accept_reject_update!(a::Assignment, swap::Swap, g, ::Strict)
    # calculate the score of the current assignment
    current_score = score(a, g)
    # perform the swap
    apply_swap!(a, swap)
    # calculate the score of the new assignment
    new_score = score(a, g)
    # if the new assignment is worse, revert the swap
    if new_score < current_score
        revert_swap!(a, swap)
    end
end
