function accept_reject_update!(history::MultivalueHistory,
                               iteration::Int,
                               proposal::Assignment,
                               current::Assignment;)
    # fill in

    push!(history, :likelihood, iteration, current.likelihood)
    return current
end
