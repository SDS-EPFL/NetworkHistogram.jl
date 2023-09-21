@testset "oracle bandwidth test" begin
    A = [0 0 1 0 1 0 1 1 0 1
        0 0 1 1 1 1 1 1 0 0
        1 1 0 1 0 0 0 0 1 0
        0 1 1 0 1 0 1 0 0 0
        1 1 0 1 0 0 1 0 0 1
        0 1 0 0 0 0 0 1 0 0
        1 1 0 1 1 0 0 1 0 1
        1 1 0 0 0 1 1 0 0 1
        0 0 1 0 0 0 0 0 0 1
        1 0 0 0 1 0 1 1 1 0]
    h = NetworkHistogram.oracle_bandwidth(A)
    rho = sum(A) / (size(A, 1) * (size(A, 1) - 1))
    h_true_nethist = 2.643731 # version 0.2.3 from nethist package
    h_clean = 3
    @test h≈h_true_nethist atol=1e-4
    @test NetworkHistogram.select_bandwidth(A) == h_clean
end
