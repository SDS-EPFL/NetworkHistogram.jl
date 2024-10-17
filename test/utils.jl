function to_default_assignment(a_specialised::Assignment{T, B}) where {T, B}
    return Assignment(a_specialised.group_size, a_specialised.node_labels)
end

to_default_assignment(a::Assignment{T, Nothing}) where {T} = a
