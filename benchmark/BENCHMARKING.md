# NetworkHistogram Performance Benchmarking Suite

This directory contains tools for measuring and tracking the performance of
NetworkHistogram's optimization algorithms.

## Overview

The benchmarking suite is designed to:

1. **Track performance improvements/regressions** over time
2. **Identify bottlenecks** in the optimization workflow
3. **Ensure optimization changes maintain correctness**
4. **Compare performance** before and after code changes

## Files

- `benchmark_optimization.jl` - Standalone benchmarking script that runs all
  benchmarks
- `test_performance_regression.jl` - Test suite that can be run with
  `Pkg.test()`
- `benchmark_results/` - Directory for storing benchmark results (created
  automatically)

## Quick Start

### Running Benchmarks

```bash
# From the repository root
julia --project=. benchmark_optimization.jl

# With custom output file
julia --project=. benchmark_optimization.jl results/my_benchmark.json

# Compare with baseline
julia --project=. benchmark_optimization.jl results/current.json results/baseline.json
```

### Running as Tests

```bash
# Run all tests including performance tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run only performance tests
julia --project=test test/test_performance_regression.jl
```

## Workflow for Performance Optimization

### 1. Establish Baseline

Before making any changes, establish a baseline:

```bash
julia --project=. benchmark_optimization.jl benchmark_results/baseline.json
```

### 2. Make Your Changes

Edit the source code to improve performance (e.g., optimize `apply_swap!`,
reduce allocations, etc.).

### 3. Run Benchmarks

```bash
julia --project=. benchmark_optimization.jl benchmark_results/after_changes.json
```

### 4. Compare Results

The script will automatically compare with `baseline.json` if it exists, or you
can manually compare:

```bash
julia --project=. benchmark_optimization.jl \
    benchmark_results/after_changes.json \
    benchmark_results/baseline.json
```

### 5. Verify Correctness

Run the full test suite to ensure your changes don't break anything:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Benchmark Categories

### Single Swap Operations

Measures the performance of a single node swap operation (apply + revert):

- **Bernoulli networks**: Binary edge weights (0/1)

  - Small: n=50, k=2
  - Medium: n=200, k=3
  - Large: n=500, k=5

- **Categorical networks**: Multi-valued edge weights
  - Small: n=50, k=2, m=3
  - Medium: n=200, k=3, m=4
  - Large: n=500, k=5, m=5

**Why it matters**: Swap operations are the core of the greedy optimization
algorithm and are called millions of times.

### Full Optimization Workflow

Measures end-to-end performance of the optimization process:

- Bernoulli: n=100, 1,000 iterations
- Categorical: n=100, 1,000 iterations

**Why it matters**: Shows real-world performance for typical use cases.

### Component Benchmarks

Measures individual components:

- **Assignment creation**: Time to create initial assignment
- **EdgeList creation**: Time to convert adjacency matrix to edge list
- **Loglikelihood computation**: Time to compute total log-likelihood
- **Get edges in groups**: Time to extract edges between two groups

**Why it matters**: Identifies which components are bottlenecks.

## Interpreting Results

### Benchmark Output

```
Benchmarking Bernoulli swap (n=50, k=2)...
  Median: 0.234 ms
```

- **Median**: The middle value (most representative of typical performance)
- **Mean**: Average value (affected by outliers)
- **Min/Max**: Best and worst case performance
- **Std**: Standard deviation (consistency of performance)

### Performance Comparison

```
✓ FASTER bernoulli_swap_n50_k2: 1.23x (23.0%)
         Current: 0.190 ms | Baseline: 0.234 ms

✗ SLOWER categorical_swap_n200_k3_m4: 0.87x (-13.0%)
         Current: 1.450 ms | Baseline: 1.260 ms

≈ SIMILAR bernoulli_optimize_n100_1k: 1.02x (2.0%)
         Current: 123.4 ms | Baseline: 125.9 ms
```

- **✓ FASTER**: >5% improvement
- **✗ SLOWER**: >5% regression
- **≈ SIMILAR**: Within ±5%

## Key Performance Hotspots

Based on the codebase analysis, these are the most critical areas for
optimization:

### 1. `apply_swap!` (swap_workspace.jl, swap_categorical.jl)

**Impact**: Called once per iteration in greedy search

**Current approach**:

- Iterates over all neighbors of swapped nodes
- Updates θ parameters and log-likelihoods incrementally
- Uses `deepcopy` for categorical distributions

**Optimization opportunities**:

- Reduce allocations in the hot path
- Optimize neighbor iteration
- Pre-allocate workspace for intermediate computations

### 2. `get_edges_in_groups` (assignment.jl)

**Impact**: Called during likelihood recomputation

**Current approach**:

- Allocates new vector for each call
- Uses `findall` to identify nodes in groups
- Iterates over all edges

**Optimization opportunities**:

- Pre-compute and cache group membership
- Use pre-allocated buffers
- Use views instead of copying data

### 3. Edge iteration (EdgeList.jl)

**Impact**: Used throughout the codebase

**Current approach**:

- Iterator-based access to edges

**Optimization opportunities**:

- Ensure type stability
- Minimize bounds checking
- Cache frequently accessed data

### 4. Log-likelihood computation

**Impact**: Called after every swap

**Current approach**:

- Recomputes for affected groups only (good!)
- Calls `logpdf` for each edge

**Optimization opportunities**:

- Batch logpdf computations
- Use SIMD operations where possible
- Cache intermediate results

## Example: Optimizing a Function

Let's say you want to optimize `apply_swap!`:

```julia
# 1. Add profiling annotations
using Profile

@profile begin
    for i in 1:1000
        apply_swap!(assignment, swap)
        revert_swap!(assignment, swap)
    end
end

Profile.print()

# 2. Identify hot spots from profiling output

# 3. Make targeted changes (e.g., reduce allocations)

# 4. Benchmark before and after
julia benchmark_optimization.jl
```

## Tips for Performance Optimization

1. **Start with profiling**: Use `@profile` to identify actual bottlenecks
2. **Benchmark incrementally**: Make one change at a time
3. **Check allocations**: Use `@btime` with `samples=1 evals=1` to see
   allocations
4. **Maintain correctness**: Always run tests after changes
5. **Consider trade-offs**: Sometimes slight speedups aren't worth added
   complexity

## Advanced Usage

### Custom Benchmarks

Add your own benchmarks to `benchmark_optimization.jl`:

```julia
function benchmark_my_function()
    # Setup
    data = create_test_data()

    # Benchmark
    b = @benchmark my_function($data) samples=100

    return Dict(
        "median_ms" => median(b.times) / 1e6,
        "mean_ms" => mean(b.times) / 1e6
    )
end
```

### Continuous Integration

To track performance over time in CI:

```yaml
# .github/workflows/benchmark.yml
- name: Run benchmarks
  run: julia --project=. benchmark_optimization.jl results/current.json

- name: Compare with main
  run: |
    git checkout main
    julia --project=. benchmark_optimization.jl results/baseline.json
    git checkout -
    julia --project=. benchmark_optimization.jl results/current.json results/baseline.json
```

## Troubleshooting

### Inconsistent Results

If you see high variance in results:

- Close other applications
- Run with `--threads=1` to avoid threading variability
- Increase the number of samples
- Let the system warm up with a few iterations first

### Out of Memory

For large benchmarks:

- Reduce the number of samples
- Run benchmarks separately instead of all at once
- Use smaller test networks

### Compilation Effects

Julia's JIT compilation can affect first-run timing:

- BenchmarkTools automatically handles warmup
- For manual timing, always run at least once before measuring

## Resources

- [BenchmarkTools.jl documentation](https://juliaci.github.io/BenchmarkTools.jl/stable/)
- [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [Profile module documentation](https://docs.julialang.org/en/v1/stdlib/Profile/)
