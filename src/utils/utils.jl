function _binary_to_index(binary_vector::Vector{Int})
    total = 1
    for i in 1:length(binary_vector)
        total += binary_vector[i] * 2^(i - 1)
    end
    return total
end

function _index_to_binary(index::Int, M)
    binary_vector = zeros(Int, M)
    index -= 1
    for i in 1:M
        binary_vector[i] = index % 2
        index = index ÷ 2
    end
    return binary_vector
end
