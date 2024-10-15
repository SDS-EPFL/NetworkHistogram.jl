abstract type AcceptRule end
struct Strict <: AcceptRule end

"""
    accept_reject_update!(proposal::Assignment, current::Assignment,
                          accept_rule::AcceptRule)


Return the updated `current` assignment based on the `accept_rule`.

# Implemented rules
- `Strict()`: Accept the proposal if it has a higher likelihood than the current assignment.
"""
accept_reject_update!

function accept_reject_update!(swap::Swap, current::Assignment, ::Strict)
end
