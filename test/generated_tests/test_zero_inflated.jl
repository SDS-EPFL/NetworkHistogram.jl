using Test
using Distributions
using Random
using NetworkHistogram: ZeroInflated, get_proba_zero

@testset "ZeroInflated Distribution Tests" begin
    @testset "continuous distribution" begin
        # Test construction
        dist = Normal(0, 1)
        zero_inflated_dist = ZeroInflated(0.5, dist)
        @test zero_inflated_dist.edge_proba == Bernoulli(0.5)
        @test zero_inflated_dist.dist == dist

        # Test pdf
        @test pdf(zero_inflated_dist, 0) ≈ 0.5 + 0.5 * pdf(dist, 0)
        @test pdf(zero_inflated_dist, 1) ≈ 0.5 * pdf(dist, 1)

        # Test get_proba_zero
        @test get_proba_zero(zero_inflated_dist) == 0.5

        # Test rand
        rng = MersenneTwister(1234)
        sample = rand(rng, zero_inflated_dist)
        @test sample == 0 || insupport(dist, sample)

        # Test logpdf
        @test logpdf(zero_inflated_dist, 0) ≈ log(0.5* (1 + pdf(dist, 0)))
        @test logpdf(zero_inflated_dist, 1) ≈ log(0.5 * pdf(dist, 1))

        # Test minimum and maximum
        @test minimum(zero_inflated_dist) == minimum(dist)
        @test maximum(zero_inflated_dist) == maximum(dist)

        # Test insupport
        @test insupport(zero_inflated_dist, 0)
        @test insupport(zero_inflated_dist, 1) == insupport(dist, 1)

        # Test cdf
        @test cdf(zero_inflated_dist, 0) ≈ 0.5 + 0.5 * cdf(dist, 0)
        @test cdf(zero_inflated_dist, 1) ≈ 0.5 + 0.5 * cdf(dist, 1)

        # Test params
        @test params(zero_inflated_dist) == (0.5, params(dist)...)

        # Test fit
        data = [0, 0, 1, 2, 3]
        fitted_dist = fit(ZeroInflated{Bernoulli, Normal}, data, 2)
        @test fitted_dist.edge_proba == Bernoulli(0.6)
        @test fitted_dist.dist isa Normal
    end

    @testset "discrete distribution" begin
        # Test construction with discrete distribution
        dist_disc = Poisson(3)
        zero_inflated_dist_disc = ZeroInflated(0.5, dist_disc)
        @test zero_inflated_dist_disc.edge_proba == Bernoulli(0.5)
        @test zero_inflated_dist_disc.dist == dist_disc

        # Test pdf with discrete distribution
        @test pdf(zero_inflated_dist_disc, 0) ≈ 0.5 + 0.5 * pdf(dist_disc, 0)
        @test pdf(zero_inflated_dist_disc, 1) ≈ 0.5 * pdf(dist_disc, 1)

        # Test get_proba_zero with discrete distribution
        @test get_proba_zero(zero_inflated_dist_disc) == 0.5

        # Test rand with discrete distribution
        rng = MersenneTwister(1234)
        sample_disc = rand(rng, zero_inflated_dist_disc)
        @test sample_disc == 0 || insupport(dist_disc, sample_disc)

        # Test logpdf with discrete distribution
        @test logpdf(zero_inflated_dist_disc, 0) ≈ log(0.5 * (1 + pdf(dist_disc, 0)))
        @test logpdf(zero_inflated_dist_disc, 1) ≈ log(0.5 * pdf(dist_disc, 1))

        # Test minimum and maximum with discrete distribution
        @test minimum(zero_inflated_dist_disc) == minimum(dist_disc)
        @test maximum(zero_inflated_dist_disc) == maximum(dist_disc)

        # Test insupport with discrete distribution
        @test insupport(zero_inflated_dist_disc, 0)
        @test insupport(zero_inflated_dist_disc, 1) == insupport(dist_disc, 1)

        # Test cdf with discrete distribution
        @test cdf(zero_inflated_dist_disc, 0) ≈ 0.5 + 0.5 * cdf(dist_disc, 0)
        @test cdf(zero_inflated_dist_disc, 1) ≈ 0.5 + 0.5 * cdf(dist_disc, 1)

        # Test params with discrete distribution
        @test params(zero_inflated_dist_disc) == (0.5, params(dist_disc)...)

        # Test fit with discrete distribution
        data_disc = [0, 0, 1, 2, 3]
        fitted_dist_disc = fit(ZeroInflated{Bernoulli, Poisson}, data_disc, 2)
        @test fitted_dist_disc.edge_proba == Bernoulli(0.6)
        @test fitted_dist_disc.dist isa Poisson
    end
end
