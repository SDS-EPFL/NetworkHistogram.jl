function accept_reject_update!(history, i, proposal, current; )
    # fill in
    
    push!(history, :likelihood, i, current.likelihood)
    return current
end
