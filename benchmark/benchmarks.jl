# using BenchmarkTools

# SUITE = BenchmarkGroup()
# for file in readdir(@__DIR__)
#     if startswith(file, "bench_") && endswith(file, ".jl")
#         SUITE[file[length("bench_") + 1:end - length(".jl")]] =
#             include(file)
#     end
# end

using BenchmarkTools, Random, Distributions, LinearAlgebra
import NetworkHistogram as NH
const SUITE = BenchmarkGroup()

function make_A(n, dist)
    A = zeros(Int, n, n)
    for j in 1:n
        for i in j:n
            if i == j
                A[i, j] = 0
            else
                A[i, j] = rand(dist)
                A[j, i] = A[i, j]
            end
        end
    end
    return A
end
# Create hierarchy of benchmarks:
SUITE["Assignment"] = BenchmarkGroup(["assignment"])

Random.seed!(123451)
stop_rule = NH.PreviousBestValue(200)
iterations = 200
swap_rule = NH.RandomNodeSwap()
accept_rule = NH.Strict()
dist = Bernoulli(0.5)

for ae in ["Bernoulli", "default"]
    if ae == "default"
        init_rule = NH.InitRule(NH.OrderedStart(), nothing)
    else
        init_rule = NH.InitRule(NH.OrderedStart(), Val{NH.BernoulliData}())
    end
    for n in [60,120,300]
        obs = NH.Observations(make_A(n,dist), dist)
        h = n ÷ 20
        a = NH.make_assignment(obs, h, init_rule)
        swap = NH.make_swap(a, (1, n))
        SUITE["Assignment"][ae]["local_search!"][n] = @benchmarkable NH.local_search!(
            $a, $obs, $swap, swap_rule = $swap_rule, accept_rule = $accept_rule)
    end
end

# tune!(SUITE);
results = run(SUITE, verbose = true)
