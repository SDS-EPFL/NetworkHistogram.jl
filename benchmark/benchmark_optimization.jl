"""
Standalone benchmarking script for NetworkHistogram optimization.

This script runs performance benchmarks and saves results to track improvements
over time. Results are saved in JSON format with timestamps.

Usage:
    julia --project=. benchmark/benchmark_optimization.jl [output_file]

If output_file is not provided, results are saved to:
    benchmark/benchmark_results/benchmark_YYYY-MM-DD_HH-MM-SS.json
"""

using Random
using StatsBase
using StaticArrays
using BenchmarkTools
using JSON3
using Dates
using NetworkHistogram

# Create output directory if it doesn't exist
const BENCHMARK_DIR = joinpath(@__DIR__, "benchmark_results")
mkpath(BENCHMARK_DIR)

# Helper functions to create test networks
function create_test_sbm_bernoulli(n_groups::Int, n_nodes::Int; seed = 42)
    Random.seed!(seed)
    d = NetworkHistogram.Bernoulli(0.5)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    for g1 in 1:n_groups
        for g2 in g1:n_groups
            p = 0.1 + 0.7 * rand()
            sbm[g1, g2] = NetworkHistogram.Bernoulli(p)
        end
    end

    base_size = n_nodes ÷ n_groups
    remainder = n_nodes % n_groups
    sizes = fill(base_size, n_groups)
    sizes[1:remainder] .+= 1
    labels = StatsBase.inverse_rle(1:n_groups, sizes)
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

function create_test_sbm_categorical(
        n_groups::Int, n_nodes::Int, n_categories::Int; seed = 42)
    Random.seed!(seed)
    ps = SVector{n_categories}(fill(1 / n_categories, n_categories))
    d = NetworkHistogram.Cat(ps)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    for g1 in 1:n_groups
        for g2 in g1:n_groups
            probs = rand(n_categories)
            probs ./= sum(probs)
            sbm[g1, g2] = NetworkHistogram.Cat(SVector{n_categories}(probs))
        end
    end

    labels = StatsBase.inverse_rle(1:n_groups, fill(n_nodes ÷ n_groups, n_groups))
    # Ensure we have exactly n_nodes by padding with last group if needed
    while length(labels) < n_nodes
        push!(labels, n_groups)
    end
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

function benchmark_single_swap(
        network_type, n_nodes, n_groups, n_categories = nothing; samples = 100)
    if network_type == :bernoulli
        A, labels, d = create_test_sbm_bernoulli(n_groups, n_nodes)
    else
        A, labels, d = create_test_sbm_categorical(n_groups, n_nodes, n_categories)
    end

    edgelist = NetworkHistogram.EdgeList(A)
    assignment = NetworkHistogram.Assignment(labels, edgelist, NetworkHistogram.Dist(d))
    swap = NetworkHistogram.make_swap(assignment, (1, n_nodes))

    b = @benchmark begin
        NetworkHistogram.apply_swap!($assignment, $swap)
        NetworkHistogram.revert_swap!($assignment, $swap)
    end setup=(NetworkHistogram.make_swap_workspace!($swap.workspace, $assignment)) samples=samples evals=1

    return Dict(
        "median_ms" => median(b.times) / 1e6,
        "mean_ms" => mean(b.times) / 1e6,
        "min_ms" => minimum(b.times) / 1e6,
        "max_ms" => maximum(b.times) / 1e6,
        "std_ms" => std(b.times) / 1e6
    )
end

function benchmark_full_optimization(
        network_type, n_nodes, n_groups, n_categories = nothing,
        max_iter = 1000; samples = 10)
    if network_type == :bernoulli
        A, labels, d = create_test_sbm_bernoulli(n_groups, n_nodes)
    else
        A, labels, d = create_test_sbm_categorical(n_groups, n_nodes, n_categories)
    end

    initial_labels = rand(1:n_groups, n_nodes)

    b = @benchmark begin
        # Create fresh params for each benchmark iteration
        params = NetworkHistogram.GreedyParams(
            $max_iter,
            NetworkHistogram.RandomNodeSwap(),
            NetworkHistogram.Strict(),
            NetworkHistogram.PreviousBestValue($max_iter ÷ 2),
            false
        )
        NetworkHistogram.nethist($A, $d, $initial_labels, params)
    end samples=samples evals=1

    return Dict(
        "median_ms" => median(b.times) / 1e6,
        "mean_ms" => mean(b.times) / 1e6,
        "min_ms" => minimum(b.times) / 1e6,
        "max_ms" => maximum(b.times) / 1e6,
        "std_ms" => std(b.times) / 1e6
    )
end

function benchmark_component(component_name, setup_fn, benchmark_fn; samples = 1000)
    setup_data = setup_fn()

    b = @benchmark $benchmark_fn($setup_data...) samples=samples

    return Dict(
        "median_us" => median(b.times) / 1e3,
        "mean_us" => mean(b.times) / 1e3,
        "min_us" => minimum(b.times) / 1e3,
        "max_us" => maximum(b.times) / 1e3,
        "std_us" => std(b.times) / 1e3
    )
end

function run_all_benchmarks()
    println("="^70)
    println("NetworkHistogram Performance Benchmarks")
    println("Started at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println("="^70)

    results = Dict(
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "julia_version" => string(VERSION),
        "benchmarks" => Dict()
    )

    # Single swap benchmarks - Bernoulli
    println("\n--- Single Swap Operations (Bernoulli) ---")
    for (n, k, s) in [(50, 2, 100), (200, 3, 100), (500, 5, 50)]
        println("Benchmarking Bernoulli swap (n=$n, k=$k)...")
        results["benchmarks"]["bernoulli_swap_n$(n)_k$(k)"] = benchmark_single_swap(
            :bernoulli, n, k; samples = s)
        r = results["benchmarks"]["bernoulli_swap_n$(n)_k$(k)"]
        println("  Median: $(round(r["median_ms"], digits=3)) ms")
    end

    # Single swap benchmarks - Categorical
    println("\n--- Single Swap Operations (Categorical) ---")
    for (n, k, m, s) in [(50, 2, 3, 100), (200, 3, 4, 100), (500, 5, 5, 50)]
        println("Benchmarking Categorical swap (n=$n, k=$k, m=$m)...")
        results["benchmarks"]["categorical_swap_n$(n)_k$(k)_m$(m)"] = benchmark_single_swap(
            :categorical, n, k, m; samples = s)
        r = results["benchmarks"]["categorical_swap_n$(n)_k$(k)_m$(m)"]
        println("  Median: $(round(r["median_ms"], digits=3)) ms")
    end

    # Full optimization benchmarks
    println("\n--- Full Optimization Workflow ---")
    println("Benchmarking Bernoulli optimization (n=100, 1k iters)...")
    results["benchmarks"]["bernoulli_optimize_n100_1k"] = benchmark_full_optimization(
        :bernoulli, 100, 3, nothing, 1000; samples = 10)
    r = results["benchmarks"]["bernoulli_optimize_n100_1k"]
    println("  Median: $(round(r["median_ms"], digits=1)) ms")

    println("Benchmarking Categorical optimization (n=100, 1k iters)...")
    results["benchmarks"]["categorical_optimize_n100_1k"] = benchmark_full_optimization(
        :categorical, 100, 3, 3, 1000; samples = 10)
    r = results["benchmarks"]["categorical_optimize_n100_1k"]
    println("  Median: $(round(r["median_ms"], digits=1)) ms")

    # Component benchmarks
    println("\n--- Component Benchmarks ---")

    println("Benchmarking Assignment creation (n=200)...")
    results["benchmarks"]["assignment_creation_n200"] = benchmark_component(
        "assignment_creation",
        () -> begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)
            (labels, edgelist, NetworkHistogram.Dist(d))
        end,
        (labels, edgelist, dist) -> NetworkHistogram.Assignment(labels, edgelist, dist);
        samples = 100
    )
    r = results["benchmarks"]["assignment_creation_n200"]
    println("  Median: $(round(r["median_us"], digits=1)) μs")

    println("Benchmarking EdgeList creation (n=200)...")
    results["benchmarks"]["edgelist_creation_n200"] = benchmark_component(
        "edgelist_creation",
        () -> begin
            A, _, _ = create_test_sbm_bernoulli(3, 200)
            (A,)
        end,
        (A,) -> NetworkHistogram.EdgeList(A);
        samples = 100
    )
    r = results["benchmarks"]["edgelist_creation_n200"]
    println("  Median: $(round(r["median_us"], digits=1)) μs")

    println("Benchmarking Loglikelihood computation (n=200)...")
    results["benchmarks"]["loglikelihood_n200"] = benchmark_component(
        "loglikelihood",
        () -> begin
            A, labels, d = create_test_sbm_bernoulli(3, 200)
            edgelist = NetworkHistogram.EdgeList(A)
            assignment = NetworkHistogram.Assignment(
                labels, edgelist, NetworkHistogram.Dist(d))
            (assignment,)
        end,
        (assignment,) -> NetworkHistogram.loglikelihood(assignment);
        samples = 1000
    )
    r = results["benchmarks"]["loglikelihood_n200"]
    println("  Median: $(round(r["median_us"], digits=2)) μs")

    return results
end

function save_results(results, output_file = nothing)
    if isnothing(output_file)
        timestamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
        output_file = joinpath(BENCHMARK_DIR, "benchmark_$timestamp.json")
    end

    open(output_file, "w") do io
        JSON3.pretty(io, results)
    end

    println("\n" * "="^70)
    println("Results saved to: $output_file")
    println("="^70)

    return output_file
end

function compare_with_baseline(results, baseline_file)
    if !isfile(baseline_file)
        println("\nBaseline file not found: $baseline_file")
        return
    end

    baseline = JSON3.read(read(baseline_file, String))

    println("\n" * "="^70)
    println("Performance Comparison vs Baseline")
    println("Baseline: $(baseline["timestamp"])")
    println("="^70)

    for (key, value) in results["benchmarks"]
        if haskey(baseline["benchmarks"], key)
            baseline_val = baseline["benchmarks"][key]

            # Determine which unit to use (ms or us)
            # Check both string and symbol keys for JSON3 compatibility
            if haskey(value, "median_ms")
                current_median = value["median_ms"]
                # Try string key first, then symbol key
                if haskey(baseline_val, "median_ms")
                    baseline_median = baseline_val["median_ms"]
                elseif haskey(baseline_val, :median_ms)
                    baseline_median = baseline_val[:median_ms]
                else
                    baseline_median = get(
                        baseline_val, "median_us", get(baseline_val, :median_us, 0)) / 1000
                end
                unit = "ms"
            else
                current_median = value["median_us"]
                # Try string key first, then symbol key
                if haskey(baseline_val, "median_us")
                    baseline_median = baseline_val["median_us"]
                elseif haskey(baseline_val, :median_us)
                    baseline_median = baseline_val[:median_us]
                else
                    baseline_median = get(
                        baseline_val, "median_ms", get(baseline_val, :median_ms, 0)) * 1000
                end
                unit = "μs"
            end

            speedup = baseline_median / current_median
            change_pct = (speedup - 1) * 100

            status = if speedup > 1.05
                "✓ FASTER"
            elseif speedup < 0.95
                "✗ SLOWER"
            else
                "≈ SIMILAR"
            end

            println("$status $key: $(round(speedup, digits=2))x ($(round(change_pct, sigdigits=3))%)")
            println("         Current: $(round(current_median, digits=2)) $unit | Baseline: $(round(baseline_median, digits=2)) $unit")
        end
    end
end

# Main execution
function main()
    output_file = length(ARGS) >= 1 ? ARGS[1] : nothing
    baseline_file = length(ARGS) >= 2 ? ARGS[2] : joinpath(BENCHMARK_DIR, "baseline.json")

    results = run_all_benchmarks()
    saved_file = save_results(results, output_file)

    if isfile(baseline_file)
        compare_with_baseline(results, baseline_file)
    else
        println("\nNo baseline found. To set this as baseline, run:")
        println("  cp $saved_file $baseline_file")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
