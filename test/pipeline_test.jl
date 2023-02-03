@testset begin "Pipeline"
    A = [
        0  0  1  0  1  0  1  1  0  1
        0  0  1  1  1  1  1  1  0  0
        1  1  0  1  0  0  0  0  1  0
        0  1  1  0  1  0  1  0  0  0
        1  1  0  1  0  0  1  0  0  1
        0  1  0  0  0  0  0  1  0  0
        1  1  0  1  1  0  0  1  0  1
        1  1  0  0  0  1  1  0  0  1
        0  0  1  0  0  0  0  0  0  1
        1  0  0  0  1  0  1  1  1  0
    ]
    graphist = NetworkHistogram.graphhist(A; h = 0.5)
end
