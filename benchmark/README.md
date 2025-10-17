# NetworkHistogram Benchmarks

This directory contains benchmarking and profiling tools for
NetworkHistogram.jl performance analysis.

## Quick Start

```bash
# Run all benchmarks (saves to benchmark_results/ with timestamp)
julia --project=. benchmark/benchmark_optimization.jl

# Profile to find bottlenecks
julia --project=. benchmark/profile_optimization.jl swap
```

## Files

| File                        | Purpose                                               |
| --------------------------- | ----------------------------------------------------- |
| `benchmark_optimization.jl` | Main benchmarking script - runs all performance tests |
| `profile_optimization.jl`   | Profile code to identify bottlenecks                  |
| `visualize_benchmarks.jl`   | Compare benchmark results over time                   |
| `run_benchmarks.jl`         | Convenience wrapper with baseline management          |
| `benchmark_results/`        | Stored benchmark results (JSON, auto-created)         |

## Usage Examples

### Running Benchmarks

```bash
# Basic usage
julia --project=. benchmark/benchmark_optimization.jl

# Save to specific file
julia --project=. benchmark/benchmark_optimization.jl my_results.json

# Compare with baseline (auto-detects baseline.json)
julia --project=. benchmark/benchmark_optimization.jl
```

### Profiling

```bash
# Profile swap operations
julia --project=. benchmark/profile_optimization.jl swap

# Profile full optimization
julia --project=. benchmark/profile_optimization.jl optimize

# Profile individual components
julia --project=. benchmark/profile_optimization.jl components
```

### Baseline Management

```bash
# Set current run as baseline
cp benchmark/benchmark_results/benchmark_2025-10-15T22-51-21.json \
   benchmark/benchmark_results/baseline.json
```

## What Gets Benchmarked

### Single Swap Operations

Tests the core swap operation (apply + revert):

- **Bernoulli networks**: Binary edges (0/1)
  - Small: n=50, k=2
  - Medium: n=200, k=3
  - Large: n=500, k=5
- **Categorical networks**: Multi-valued edges
  - Small: n=50, k=2, m=3
  - Medium: n=200, k=3, m=4
  - Large: n=500, k=5, m=5

### Full Optimization

End-to-end optimization performance:

- Bernoulli: n=100, 1,000 iterations
- Categorical: n=100, 1,000 iterations

### Component Benchmarks

Individual function performance:

- Assignment creation
- EdgeList creation
- Log-likelihood computation
- Edge extraction between groups

## Current Performance (October 2025)

**After Phase 1 Optimizations:**

| Operation                    | Time    | vs Baseline |
| ---------------------------- | ------- | ----------- |
| Bernoulli swap (n=500)       | 2.56 ms | 6.1x faster |
| Bernoulli swap (n=200)       | 0.54 ms | 3.2x faster |
| Categorical swap (n=200)     | 0.04 ms | 1.4x faster |
| Bernoulli optimize (n=100)   | 92 ms   | -           |
| Categorical optimize (n=100) | 15 ms   | -           |

## Optimization Workflow

1. **Establish baseline**: Run benchmarks before changes
2. **Profile**: Use `profile_optimization.jl` to find hotspots
3. **Optimize**: Edit source code (usually `src/optimization/`)
4. **Benchmark**: Run benchmarks again
5. **Verify**: Run tests to ensure correctness
6. **Repeat**: Continue until satisfied

## Key Hotspots

Focus optimization efforts on:

1. **`apply_swap!`** - Called millions of times (biggest impact)
2. **`get_edges_in_groups`** - Called during likelihood updates
3. **Edge iteration** - Used throughout, cumulative effect

See `dev/CODE_REVIEW_2025-10-15.md` for detailed analysis.

## Output Format

```
--- Single Swap Operations (Bernoulli) ---
Benchmarking Bernoulli swap (n=200, k=3)...
  Median: 0.538 ms

======================================================================
Performance Comparison vs Baseline
======================================================================
✓ FASTER bernoulli_swap_n200_k3: 3.23x (223.0%)
         Current: 0.54 ms | Baseline: 1.74 ms
```

- **✓ FASTER**: >5% improvement
- **✗ SLOWER**: >5% regression
- **≈ SIMILAR**: Within ±5%

## Tips

- **Close other apps** for consistent results
- **Run multiple times** to warm up JIT compiler (BenchmarkTools handles this)
- **Check allocations** with `@btime ... samples=1 evals=1`
- **Profile first** before optimizing
- **Test after** every optimization

## Documentation

- Full details: See `PERFORMANCE.md` (root directory)
- Code review: See `dev/CODE_REVIEW_2025-10-15.md`
