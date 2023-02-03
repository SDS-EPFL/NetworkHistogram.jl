function accept_reject_update!(history::MVHistory, iteration::Int, proposal::Assignment,
    current::Assignment; accept_rule)
# fill in

push!(history, :likelihood, iteration, current.likelihood)
return current
end
