# NetworkHistogram Benchmarks

This directory contains benchmarking and profiling tools for
NetworkHistogram.jl performance analysis.

## Files

### `benchmark_optimization.jl`

Main benchmarking script that runs comprehensive performance tests.

**Usage:**

```bash
julia --project=. benchmark/benchmark_optimization.jl [output_file]
```

**Features:**

- Single swap operations (Bernoulli and Categorical networks)
- Full optimization workflows
- Component benchmarks (Assignment, EdgeList, loglikelihood)
- Automatic comparison with baseline
- Results saved as timestamped JSON files

**Example:**

```bash
# Run benchmarks and save to default location
julia --project=. benchmark/benchmark_optimization.jl

# Save to specific file
julia --project=. benchmark/benchmark_optimization.jl my_results.json

# Compare with custom baseline
julia --project=. benchmark/benchmark_optimization.jl current.json baseline.json
```

### `visualize_benchmarks.jl`

Compare and visualize benchmark results over time.

**Usage:**

```bash
julia --project=. benchmark/visualize_benchmarks.jl [files...]
```

**Example:**

```bash
# Compare two benchmark runs
julia --project=. benchmark/visualize_benchmarks.jl \
    benchmark/benchmark_results/baseline.json \
    benchmark/benchmark_results/benchmark_2025-10-15T22-51-21.json

# Use all files in benchmark_results/
julia --project=. benchmark/visualize_benchmarks.jl --all
```

### `profile_optimization.jl`

Profile code to identify performance bottlenecks.

**Usage:**

```bash
julia --project=. benchmark/profile_optimization.jl [scenario]
```

**Scenarios:**

- `swap` - Profile single swap operations
- `optimize` - Profile full optimization run
- `components` - Profile individual components

**Example:**

```bash
julia --project=. benchmark/profile_optimization.jl swap
```

## Benchmark Results

Results are stored in `benchmark/benchmark_results/` as JSON files with
timestamps.

### Setting a Baseline

To set a benchmark run as the baseline for future comparisons:

```bash
cp benchmark/benchmark_results/benchmark_2025-10-15T22-51-21.json \
   benchmark/benchmark_results/baseline.json
```

## Performance Metrics

The benchmarks track:

- **Median time**: Most representative performance measurement
- **Mean time**: Average across all samples
- **Min/Max time**: Best and worst case performance
- **Standard deviation**: Performance consistency

### Current Performance (October 2025)

**Bernoulli Swap Operations:**

- n=50, k=2: ~0.038 ms
- n=200, k=3: ~0.52 ms
- n=500, k=5: ~2.4 ms (6.3x faster than pre-optimization baseline)

**Categorical Swap Operations:**

- n=50, k=2, m=3: ~0.009 ms
- n=200, k=3, m=4: ~0.04 ms
- n=500, k=5, m=5: ~0.10 ms

**Full Optimization:**

- Bernoulli (n=100, 1k iters): ~90 ms
- Categorical (n=100, 1k iters): ~20 ms

## Optimization History

Major optimizations implemented:

1. Eliminated `deepcopy` in swap operations (replaced with in-place
   `copy_symarray!`)
2. Optimized `get_edges_in_groups` (removed `findall`, added direct iteration)
3. Added `@inbounds` annotations to hot paths
4. Pre-sized Set allocations to avoid resizing

**Result:** 6.3x speedup for large Bernoulli networks while maintaining full
correctness.
