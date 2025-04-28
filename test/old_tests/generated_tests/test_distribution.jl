using NetworkHistogram: ZeroInflated, DiscretizedDistribution,
                        ZeroInflatedCategorical,
                        ncategories, Discretizer, encode, decode, binwidth,
                        RegularDiscretizer,
                        CategoryDiscretizer, HybridDiscretizer,
                        DiscretizerZeroToZero, nlabels
using Distributions
using Test

@testset "ZeroInflated" begin
    dist = ZeroInflated(0.3, truncated(Normal(0, 1), -3, 3))
    @test pdf(dist, 0) ≈ 0.3 + 0.7 * pdf(truncated(Normal(0, 1), -3, 3), 0)
    @test pdf(dist, 1) ≈ 0.7 * pdf(truncated(Normal(0, 1), -3, 3), 1)
    @test cdf(dist, 0) ≈ 0.3 + 0.7 * cdf(truncated(Normal(0, 1), -3, 3), 0)
    @test cdf(dist, 1) ≈ 0.3 + 0.7 * cdf(truncated(Normal(0, 1), -3, 3), 1)
end

@testset "DiscretizedDistribution" begin
    dist = DiscretizedDistribution(truncated(Normal(0, 1), -3, 3), 10)
    @test ncategories(dist) == 10
    @test pdf(dist, 0) >= 0
    @test cdf(dist, 0) >= 0
end

@testset "ZeroInflatedCategorical" begin
    dist = ZeroInflatedCategorical(0.3, Categorical([0.2, 0.3, 0.5]))
    @test pdf(dist, 0) ≈ 0.3
    @test pdf(dist, 1) ≈ 0.7 * 0.2
    @test cdf(dist, 0) ≈ 0.3
    @test cdf(dist, 1) ≈ 0.3 + 0.7 * 0.2
end

@testset "ZeroInflatedDiscretizedDistribution" begin
    dist = ZeroInflated(0.3, truncated(Normal(0, 1), -3, 3))
    disc_dist = DiscretizedDistribution(dist, 10)
    @test ncategories(disc_dist) == 10
    @test pdf(disc_dist, 0) >= 0
    @test cdf(disc_dist, 0) >= 0
end

@testset "DiscretizedZeroInflatedCategorical" begin
    dist = ZeroInflatedCategorical(0.3, Categorical([0.2, 0.3, 0.5]))
    disc_dist = DiscretizedDistribution(dist, 10)
    @test ncategories(disc_dist) == 10
    @test pdf(disc_dist, 0) >= 0
    @test cdf(disc_dist, 0) >= 0
end

@testset "Discretizer" begin
    using Distributions
    disc = RegularDiscretizer(10, 0.0, 1.0)
    @test encode(disc, 0.05) == 1
    @test decode(disc, 1) == (0.0, 0.1)
    @test binwidth(disc) == 0.1
    @test nlabels(disc) == 10
end

@testset "CategoryDiscretizer" begin
    cat_to_bin = Dict("a" => 1, "b" => 2, "c" => 3)
    bin_to_cat = Dict(1 => "a", 2 => "b", 3 => "c")
    disc = CategoryDiscretizer(cat_to_bin, bin_to_cat)
    @test encode(disc, "a") == 1
    @test decode(disc, 1) == "a"
    @test nlabels(disc) == 3
end

@testset "HybridDiscretizer" begin
    atoms = [0.0, 1.0]
    disc = HybridDiscretizer(10, -1.0, 1.0, atoms)
    @test encode(disc, 0.0) == 11
    @test encode(disc, 0.5) == 8
    @test decode(disc, 11) == 0.0
    @test all(isapprox.(decode(disc, 8), (0.4, 0.6); atol = 1e-2))
    @test nlabels(disc) == 12
end

@testset "DiscretizerZeroToZero" begin
    disc = DiscretizerZeroToZero(10, -1.0, 1.0)
    @test encode(disc, 0.0) == 0
    @test encode(disc, 0.5) == 8
    @test decode(disc, 0) == 0.0
    @test all(isapprox.(decode(disc, 8), (0.4, 0.6); atol = 1e-2))
    @test nlabels(disc) == 11
end
