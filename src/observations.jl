# getters for observations

function get_obs(g::SimpleGraph, x::Tuple)
    return get_obs(g, x[1], x[2])
end

function get_obs(g::SimpleGraph, i::Int, j::Int)
    return convert(Bool, has_edge(g, i, j))
end
