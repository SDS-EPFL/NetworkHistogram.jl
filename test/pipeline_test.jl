@testset begin "Pipeline"
    n = 10
    A = collect(Symmetric([i == j ? 0 : rand(0:1) for i in 1:n, j in 1:n]))
    graphist = NetworkHistogram.graphhist(A; h = 0.5)
end
