using BenchmarkTools, NetworkHistogram
const SUITE = BenchmarkGroup()

# Create hierarchy of benchmarks:
SUITE["eval"] = BenchmarkGroup()

options = Options(; binary_operators = [+, -, *], unary_operators = [cos])


for n in [10, 20]
    SUITE["eval_tree_array"][n] = @benchmarkable(eval_tree_array($tree, X, $options),
        evals=10,
        samples=1000,
        setup=(X = randn(Float32, 2, $n)))
end
