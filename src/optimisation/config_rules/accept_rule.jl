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

function slow_swap!(a, swap)
    swap_node_labels!(a, swap.index1, swap.index2)
    _count_occurences!(a.additional_data, a.node_labels)
    update_ll!(a)
end

function accept_reject_update!(a::Assignment, swap::Swap, g, ::Strict)
    # calculate the score of the current assignment
    current_score = score(a, g)
    # perform the swap
    #a_star = deepcopy(a)
    #swap_star = deepcopy(swap)
    apply_swap!(a, swap)
    #slow_swap!(a_star, swap_star)
    # if any(a_star.additional_data.realized .!= a.additional_data.realized)
    #     println("The slow and fast swap functions do not agree on realized")
    #     if any(a_star.additional_data.counts .!= a.additional_data.counts)
    #         println("The slow and fast swap functions do not agree on counts")
    #     end
    #     println("a_star.additional_data.realized: ", a_star.additional_data.realized)
    #     println("a.additional_data.realized: ", a.additional_data.realized)
    #     println("a_star.additional_data.counts: ", a_star.additional_data.counts)
    #     println("a.additional_data.counts: ", a.additional_data.counts)
    #     error("The slow and fast swap functions do not agree after swapping labels ", swap.index1, " and ", swap.index2)
    # end

    # if the new assignment is worse, revert the swap
    if score(a, g) <= current_score
        revert_swap!(a, swap)
    end
end
