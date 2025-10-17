"""
Visualize benchmark results over time.

This script reads multiple benchmark result files and creates a simple
comparison table or plot showing performance trends.

Usage:
    julia --project=. benchmark/visualize_benchmarks.jl [files...]
    julia --project=. benchmark/visualize_benchmarks.jl --all  # Use all files in benchmark_results/

Example:
    julia benchmark/visualize_benchmarks.jl \\
        benchmark/benchmark_results/baseline.json \\
        benchmark/benchmark_results/current_2024-10-15.json
"""

using JSON3
using Dates
using Printf
using Statistics
using PrettyTables

function load_benchmark(filepath)
    if !isfile(filepath)
        @warn "File not found: $filepath"
        return nothing
    end

    data = JSON3.read(read(filepath, String))
    return data
end

function extract_key_metrics(benchmark_data)
    metrics = Dict{String, Float64}()

    for (name, values) in benchmark_data["benchmarks"]
        name_str = string(name)  # Convert Symbol to String
        if haskey(values, "median_ms")
            metrics[name_str] = values["median_ms"]
        elseif haskey(values, "median_us")
            metrics[name_str] = values["median_us"] / 1000.0  # Convert to ms
        end
    end

    return metrics
end

function compare_multiple(files)
    if isempty(files)
        println("No files provided")
        return
    end

    # Load all benchmarks
    benchmarks = []
    for file in files
        data = load_benchmark(file)
        if !isnothing(data)
            push!(benchmarks,
                (
                    file = basename(file),
                    timestamp = data["timestamp"],
                    metrics = extract_key_metrics(data)
                ))
        end
    end

    if isempty(benchmarks)
        println("No valid benchmark files found")
        return
    end

    # Sort by timestamp
    sort!(benchmarks, by = b -> b.timestamp)

    # Get all metric names
    all_metrics = Set{String}()
    for b in benchmarks
        union!(all_metrics, keys(b.metrics))
    end
    all_metrics = sort(collect(all_metrics))

    # Print header
    println("\n" * "="^100)
    println("Benchmark Comparison Across Versions")
    println("="^100 * "\n")

    # Prepare data for table
    table_data = []
    baseline_vals = Dict{String, Float64}()

    for metric in all_metrics
        # Skip if metric has no values
        values = [haskey(b.metrics, metric) ? b.metrics[metric] : NaN for b in benchmarks]
        if all(isnan, values)
            continue
        end

        row = Any[metric]
        baseline_val = values[1]
        baseline_vals[metric] = baseline_val

        for (i, val) in enumerate(values)
            if isnan(val)
                push!(row, "N/A")
            else
                push!(row, round(val, digits = 2))
            end
        end

        push!(table_data, row)
    end

    if isempty(table_data)
        println("No metrics to display")
        return
    end

    # Create headers
    headers = ["Benchmark"]
    for b in benchmarks
        short_name = length(b.file) > 16 ? b.file[1:13] * "..." : b.file
        push!(headers, short_name)
    end

    # Create subheaders with timestamps
    subheaders = [""]
    for b in benchmarks
        short_ts = length(b.timestamp) > 16 ? b.timestamp[1:16] : b.timestamp
        push!(subheaders, short_ts)
    end

    # Create highlighters for improvements and regressions
    # We'll color entire rows based on whether the value improved or regressed vs baseline
    hl_improvement = TextHighlighter(
        (data, i, j) -> begin
            # Check if current value (in any column after baseline) shows improvement
            baseline_idx = 2  # First value column
            baseline_val = data[i, baseline_idx]

            if j > 2 && baseline_val isa Number && baseline_val > 0
                current_val = data[i, j]
                if current_val isa Number
                    speedup = baseline_val / current_val
                    return speedup > 1.05  # >5% improvement
                end
            end
            return false
        end,
        crayon"green"
    )

    hl_regression = TextHighlighter(
        (data, i, j) -> begin
            # Check if current value (in any column after baseline) shows regression
            baseline_idx = 2  # First value column
            baseline_val = data[i, baseline_idx]

            if j > 2 && baseline_val isa Number && baseline_val > 0
                current_val = data[i, j]
                if current_val isa Number
                    speedup = baseline_val / current_val
                    return speedup < 0.95  # >5% regression
                end
            end
            return false
        end,
        crayon"red"
    )    # Convert table_data to matrix
    data_matrix = permutedims(hcat([vcat(row...) for row in table_data]...))

    # Print table
    pretty_table(
        data_matrix;
        column_labels = headers,
        highlighters = [hl_improvement, hl_regression],
        alignment = vcat(:l, fill(:r, length(benchmarks))),
        table_format = TextTableFormat(borders = text_table_borders__unicode_rounded)
    )

    println("\nLegend: Green = >5% faster, Red = >5% slower (compared to first column)")
    println("All values in milliseconds (ms)")
    println()

    # Calculate aggregate statistics
    if length(benchmarks) >= 2
        println("Overall Summary:")
        println("-" * "="^50)

        baseline = benchmarks[1]
        for i in 2:length(benchmarks)
            current = benchmarks[i]

            speedups = Float64[]
            for metric in all_metrics
                if haskey(baseline.metrics, metric) && haskey(current.metrics, metric)
                    base_val = baseline.metrics[metric]
                    curr_val = current.metrics[metric]
                    if base_val > 0 && curr_val > 0
                        push!(speedups, base_val / curr_val)
                    end
                end
            end

            if !isempty(speedups)
                median_speedup = median(speedups)
                geomean_speedup = exp(mean(log.(speedups)))
                faster_count = count(s -> s > 1.05, speedups)
                slower_count = count(s -> s < 0.95, speedups)
                similar_count = length(speedups) - faster_count - slower_count

                println("\n$(current.file) vs $(baseline.file):")
                println("  Geometric mean speedup: $(round(geomean_speedup, digits=2))x")
                println("  Median speedup: $(round(median_speedup, digits=2))x")
                println("  Benchmarks: $faster_count faster, $slower_count slower, $similar_count similar")
            end
        end
    end
end

function main()
    if length(ARGS) == 0 || ARGS[1] in ["-h", "--help", "help"]
        println("""
        Visualize NetworkHistogram Benchmark Results
        =============================================

        Usage: julia dev/visualize_benchmarks.jl [options] [files...]

        Options:
          --all         Compare all files in benchmark_results/
          -h, --help    Show this help

        Examples:
          # Compare specific files
          julia dev/visualize_benchmarks.jl \\
              benchmark_results/baseline.json \\
              benchmark_results/current.json

          # Compare all available benchmarks
          julia dev/visualize_benchmarks.jl --all
        """)
        return
    end

    files = if ARGS[1] == "--all"
        results_dir = joinpath("dev", "benchmark_results")
        if !isdir(results_dir)
            println("Error: benchmark_results directory not found")
            return
        end

        all_files = filter(f -> endswith(f, ".json"), readdir(results_dir))
        sort!([joinpath(results_dir, f) for f in all_files])
    else
        ARGS
    end

    compare_multiple(files)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
