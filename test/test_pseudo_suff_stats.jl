using Test
using NetworkHistogram
using StaticArrays
using Distributions
import NetworkHistogram as NH

function _one_hot_vector(sample::Int, num_categories::Int)
    v = zeros(Int, num_categories)
    v[sample] = 1
    return v
end

@testset "Bernoulli" begin
    @testset "score" begin
        ss = NH.BernoulliSuffStats()
        samples = [true, false, true, true, false, true, false, false, true, true]
        for s in samples
            ss = NH.add_sample(ss, s)
        end
        d = fit_mle(Bernoulli, samples)
        @test NH.score(ss) ≈ -sum(map(Base.Fix1(logpdf, d), samples))
        @test NH.to_params(ss) == d.p
    end
end

@testset "Categorical" begin
    @testset "score" begin
        ss = NH.CategoricalSuffStats(3)
        samples = [1, 2, 1, 2, 3, 1, 2, 3, 1, 2]
        s_vec = _one_hot_vector.(samples, 3)
        for s in samples
            ss = NH.add_sample(ss, s)
        end
        d = fit_mle(Categorical, samples)
        p = probs(d)
        loss = 0.0
        for s in s_vec
            loss += sum(abs2, s - p)
        end
        @test NH.score(ss) ≈ loss
        @test NH.to_params(ss) == p
    end
end
