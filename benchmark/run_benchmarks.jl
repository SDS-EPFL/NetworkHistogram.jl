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

function ensure_dependencies()
    println("Checking dependencies...")

    # Check if BenchmarkTools and JSON3 are available
    try
        @eval using BenchmarkTools
        @eval using JSON3
    catch
        println("Installing required dependencies...")
        Pkg.activate("test")
        Pkg.add(["BenchmarkTools", "JSON3"])
        Pkg.activate(".")
    end

    println("Dependencies OK ✓")
end

function run_baseline()
    ensure_dependencies()

    baseline_file = joinpath("benchmark", "benchmark_results", "baseline.json")

    if isfile(baseline_file)
        print("Baseline already exists. Overwrite? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            println("Aborted.")
            return
        end
    end

    println("\nRunning baseline benchmarks...")
    println("This may take several minutes...\n")

    run(`julia --project=. benchmark/benchmark_optimization.jl $baseline_file`)

    println("\n✓ Baseline established at: $baseline_file")
    println("\nNext steps:")
    println("  1. Make your performance improvements")
    println("  2. Run: julia run_benchmarks.jl current")
    println("  3. Review the performance comparison")
end

function run_current()
    ensure_dependencies()

    baseline_file = joinpath("benchmark", "benchmark_results", "baseline.json")

    if !isfile(baseline_file)
        println("⚠ Warning: No baseline found!")
        println("Consider running: julia run_benchmarks.jl baseline")
        println("\nContinuing anyway...\n")
    end

    timestamp = Dates.format(Dates.now(), "yyyy-mm-ddTHH-MM-SS")
    current_file = joinpath("benchmark", "benchmark_results", "current_$timestamp.json")

    println("Running current benchmarks...")
    println("This may take several minutes...\n")

    if isfile(baseline_file)
        run(`julia --project=. benchmark/benchmark_optimization.jl $current_file $baseline_file`)
    else
        run(`julia --project=. benchmark/benchmark_optimization.jl $current_file`)
    end

    println("\n✓ Results saved to: $current_file")
end

function compare_benchmarks(file1, file2)
    ensure_dependencies()

    if !isfile(file1)
        println("Error: File not found: $file1")
        return
    end

    if !isfile(file2)
        println("Error: File not found: $file2")
        return
    end

    println("Comparing benchmarks...")
    println("  Baseline: $file2")
    println("  Current:  $file1\n")

    # Re-run comparison
    run(`julia --project=. benchmark/benchmark_optimization.jl $file1 $file2`)
end

function clean_results()
    results_dir = joinpath("benchmark", "benchmark_results")

    if !isdir(results_dir)
        println("No results directory found.")
        return
    end

    files = filter(f -> endswith(f, ".json") && f != "baseline.json", readdir(results_dir))

    if isempty(files)
        println("No benchmark results to clean.")
        return
    end

    println("Found $(length(files)) benchmark result file(s):")
    for f in files
        println("  - $f")
    end

    print("\nDelete these files? (y/N): ")
    response = readline()

    if lowercase(strip(response)) == "y"
        for f in files
            rm(joinpath(results_dir, f))
        end
        println("✓ Cleaned $(length(files)) file(s)")
    else
        println("Aborted.")
    end
end

function run_tests()
    ensure_dependencies()

    println("Running full test suite...")
    Pkg.test()
end

# Main execution
function main()
    if length(ARGS) == 0 || ARGS[1] == "help" || ARGS[1] == "-h" || ARGS[1] == "--help"
        print_help()
        return
    end

    command = ARGS[1]

    if command == "baseline"
        run_baseline()
    elseif command == "current"
        run_current()
    elseif command == "compare"
        if length(ARGS) < 3
            println("Error: compare requires two file arguments")
            println("Usage: julia run_benchmarks.jl compare FILE1 FILE2")
            return
        end
        compare_benchmarks(ARGS[2], ARGS[3])
    elseif command == "clean"
        clean_results()
    elseif command == "test"
        run_tests()
    else
        println("Error: Unknown command '$command'")
        print_help()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
