function _index_to_binary(index::Int, size::Int)
    result = falses(size)
    _index_to_binary!(result, index)
    return result
end

function _index_to_binary!(binary_vector, index)
    digits!(binary_vector, index - 1; base = 2)
end

function _binary_to_index(binary_vector)
    index = 1
    @inbounds for (i, s) in enumerate(binary_vector)
        index += s * 2^(i - 1)
    end
    return index
end
