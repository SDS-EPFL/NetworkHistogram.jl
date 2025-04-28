@testset "Distribution tests" begin
    import NetworkHistogram as NH
    d1 = NH.Bernoulli(0.5)
    d2 = NH.Bernoulli(0.7)
    my_d = NH.Dist(d1)
    d_avg = NH.add_to(my_d, d2)
    @test d_avg.counts == 2
    @test d_avg.dist.p == 0.6
    d_removed = NH.remove_from(d_avg, d2)
    @test d_removed.counts == 1
    @test d_removed.dist == d1
end
