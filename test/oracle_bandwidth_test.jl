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
    @test h≈h_true_nethist atol=1e-4
    h_clean = NetworkHistogram.sanitize_bandwidth(h, size(A, 1))
    @test h_clean == 2
end


@testset "bandwidth sanitization" begin
    # ask for 50% of the nodes in each groups
    @test NetworkHistogram.sanitize_bandwidth(0.5, 10) == 5
    # ask for 100% of the nodes in each groups
    @test NetworkHistogram.sanitize_bandwidth(1, 10) == 10
end
