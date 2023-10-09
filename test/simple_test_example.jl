"""
    make_simple_example()

Makes the simple example used in many tests.
Returns A, node_labels, group_size, assignment
"""
function make_simple_example()
    A = [0 1 1 1 0 0 1 0
        1 0 1 1 0 0 0 0
        1 1 0 0 0 0 0 0
        1 1 0 0 0 0 0 1
        0 0 0 0 0 1 1 1
        0 0 0 0 1 0 1 1
        1 0 0 0 1 1 0 0
        0 0 0 1 1 1 0 0]
    node_labels = [1, 1, 1, 1, 2, 2, 2, 2]
    group_size = NetworkHistogram.GroupSize(8, 4)
    assignment = NetworkHistogram.Assignment(A, node_labels, group_size)
    return A, node_labels, group_size, assignment
end


function make_multivariate_example()
    A = [0 1 1 1 0 0 1 0
        1 0 1 1 0 0 0 0
        1 1 0 0 0 0 0 0
        1 1 0 0 0 0 0 1
        0 0 0 0 0 1 1 1
        0 0 0 0 1 0 1 1
        1 0 0 0 1 1 0 0
        0 0 0 1 1 1 0 0]
    A = cat(A, abs.(A .- 1), dims = 3)
    for i in size(A,1)
        A[i,i,:] .= 0
    end
    node_labels = [1, 1, 1, 1, 2, 2, 2, 2]
    group_size = NetworkHistogram.GroupSize(8, 4)
    assignment = NetworkHistogram.Assignment(A, node_labels, group_size)
    return A, node_labels, group_size, assignment
end



function make_er(n,p)
    A = zeros(Int64, n, n)
    for i in 1:n
        for j in i+1:n
            if rand() < p
                A[i,j] = 1
                A[j,i] = 1
            end
        end
    end
    return A
end
