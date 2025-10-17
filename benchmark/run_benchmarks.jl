#!/usr/bin/env julia

"""
Quick-start script for running NetworkHistogram benchmarks.

Usage:
    ./run_benchmarks.jl [command] [options]

Commands:
    baseline    - Establish a baseline benchmark
    current     - Run current benchmarks (compares with baseline if available)
    compare     - Compare two benchmark files
    clean       - Remove all benchmark results
    help        - Show this help message

Examples:
    ./run_benchmarks.jl baseline
    ./run_benchmarks.jl current
    ./run_benchmarks.jl compare results1.json results2.json
"""

using Pkg
using Dates

# Ensure we're in the right directory
cd(dirname(@__DIR__))

function print_help()
    println("""
    NetworkHistogram Benchmark Runner
    ==================================

    Usage: julia run_benchmarks.jl [command] [options]

    Commands:
      baseline          Create a baseline benchmark
      current           Run current benchmarks
      compare FILE1 FILE2   Compare two benchmark files
      clean             Remove all benchmark results
      help              Show this message

    Examples:
      julia run_benchmarks.jl baseline
      julia run_benchmarks.jl current
      julia run_benchmarks.jl compare results/v1.json results/v2.json
    """)
end

function ensure_dependencies(tries = 0)
    tries == 0 && @info "Checking dependencies..."
    # Check if benchmark dependencies are installed
    try
        @eval using StaticArrays
        @eval using BenchmarkTools
        @eval using JSON3
        @eval using PrettyTables
        @eval using LoggingExtras

        @info "Dependencies OK ✓"
    catch
        if tries >= 2
            error("Failed to install dependencies after multiple attempts.")
        elseif tries == 1
            @info "Trying to instantiate project..."
            Pkg.instantiate()
            ensure_dependencies(tries + 1)
        else
            @info "Activating benchmark project"
            Pkg.activate("benchmark")
            ensure_dependencies(tries + 1)
        end
    end
end

function run_baseline()
    baseline_file = joinpath("benchmark", "benchmark_results", "baseline.json")

    if isfile(baseline_file)
        printstyled(
            "Baseline already exists. Overwrite? (y/N): ", color = :light_yellow, blink = true)
        response = readline()
        if lowercase(strip(response)) != "y"
            @info "Aborted."
            return
        end
    end

    @info "Running baseline benchmarks... \nThis may take several minutes...\n"

    run(`julia --project=benchmark benchmark/benchmark_optimization.jl $baseline_file`)

    @info "\n✓ Baseline established at: $baseline_file" *
          "\n     Next steps: " *
          "\n  1.  Make your performance improvements" *
          "\n  2. Run: julia run_benchmarks.jl current" *
          "\n  3. Review the performance comparison"
end

function run_current()
    baseline_file = joinpath("benchmark", "benchmark_results", "baseline.json")

    if !isfile(baseline_file)
        @warn "⚠ Warning: No baseline found! \n Consider running: julia run_benchmarks.jl baseline"
        @info "\nContinuing anyway...\n"
    end

    timestamp = Dates.format(Dates.now(), "yyyy-mm-ddTHH-MM-SS")
    current_file = joinpath("benchmark", "benchmark_results", "current_$timestamp.json")

    @info "Running current benchmarks... \n This may take several minutes...\n"

    if isfile(baseline_file)
        run(`julia --project=benchmark benchmark/benchmark_optimization.jl $current_file $baseline_file`)
    else
        run(`julia --project=benchmark benchmark/benchmark_optimization.jl $current_file`)
    end

    @info "✓ Results saved to: $current_file"
end

function compare_benchmarks(file1, file2)
    if !isfile(file1)
        @error "Error: File not found: $file1"
        return
    end

    if !isfile(file2)
        @error "Error: File not found: $file2"
        return
    end

    @info "Comparing benchmarks... \n  Baseline: $file2 \n  Current:  $file1\n"

    # Re-run comparison
    run(`julia --project=benchmark benchmark/benchmark_optimization.jl $file1 $file2`)
end

function clean_results()
    results_dir = joinpath("benchmark", "benchmark_results")

    if !isdir(results_dir)
        @info "No results directory found."
        return
    end

    files = filter(f -> endswith(f, ".json") && f != "baseline.json", readdir(results_dir))

    if isempty(files)
        @info "No benchmark results to clean."
        return
    end

    @info "Found $(length(files)) benchmark result file(s):"
    for f in files
        @info "  - $f"
    end

    printstyled("\nDelete these files? (y/N): ", color = :light_yellow, blink = true)
    response = readline()

    if lowercase(strip(response)) == "y"
        for f in files
            rm(joinpath(results_dir, f))
        end
        @info "✓ Cleaned $(length(files)) file(s)"
    else
        @info "Aborted."
    end
end

# Main execution
function main()
    if length(ARGS) == 0 || ARGS[1] == "help" || ARGS[1] == "-h" || ARGS[1] == "--help"
        print_help()
        return
    end

    ensure_dependencies()

    command = ARGS[1]

    if command == "baseline"
        run_baseline()
    elseif command == "current"
        run_current()
    elseif command == "compare"
        if length(ARGS) < 3
            @error "Error: compare requires two file arguments \n Usage: julia run_benchmarks.jl compare FILE1 FILE2"
            return
        end
        compare_benchmarks(ARGS[2], ARGS[3])
    elseif command == "clean"
        clean_results()
    elseif command == "test"
        run_tests()
    else
        @info "Error: Unknown command '$command'"
        print_help()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
