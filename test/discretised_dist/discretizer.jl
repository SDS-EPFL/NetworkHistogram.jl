using NetworkHistogram

@testset "discretizer" begin
    using StaticArrays
    reg_disc = NetworkHistogram.RegularDiscretizer(
        10, 0.0, 1.0, MVector{10}(1:10), 1 / 10)
    cat_disc = NetworkHistogram.CategoryDiscretizer(
        Dict([0.0 => 11]), Dict([11 => 0.0]))
    hybrid_disc = NetworkHistogram.HybridDiscretizer(
        reg_disc, cat_disc)

    @test NetworkHistogram.encode(reg_disc, 0.0) == 1
    @test NetworkHistogram.encode(cat_disc, 0.0) == 11
    @test NetworkHistogram.encode(hybrid_disc, 0.0) == 11
    @test NetworkHistogram.decode(hybrid_disc, 11) == 0.0
    @test all(NetworkHistogram.encode(reg_disc, 0.001:0.001:1.0) .==
              NetworkHistogram.encode(hybrid_disc, 0.001:0.001:1.0))
    @test all(NetworkHistogram.decode(hybrid_disc, 1:10) .==
              NetworkHistogram.decode(reg_disc, 1:10))
end
