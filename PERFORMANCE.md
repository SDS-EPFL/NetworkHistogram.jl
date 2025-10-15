# Performance Optimization Guide for NetworkHistogram

This repository now includes a comprehensive performance regression test suite
to help improve the optimization speed of NetworkHistogram algorithms.

## 🎯 Quick Start

### 1. Establish a Baseline

Before making any changes:

```bash
julia dev/run_benchmarks.jl baseline
```

### 2. Make Your Changes

Edit the optimization code (e.g., in `src/optimization/`)

### 3. Test Performance

```bash
julia dev/run_benchmarks.jl current
```

This will automatically compare against your baseline and show:

- Which operations got faster/slower
- By how much (speedup factor and percentage)
- Detailed timing statistics

### 4. Verify Correctness

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## 📊 What Gets Benchmarked

### Core Operations

- **Single swap operations** (Bernoulli & Categorical networks)
  - Small networks: n=50 nodes
  - Medium networks: n=200 nodes
  - Large networks: n=500 nodes

### Full Workflows

- Complete optimization runs (1,000 iterations)
- End-to-end performance measurement

### Components

- Assignment creation
- EdgeList creation
- Log-likelihood computation
- Edge extraction

## 🔍 Key Hotspots for Optimization

Based on the workflow in `test_decorated_paper.jl`, these are the critical
bottlenecks:

### 1. `apply_swap!` Function

**Location**: `src/optimization/swap_workspace.jl`, `swap_categorical.jl`

**Why it matters**: Called millions of times during optimization (once per
iteration)

**Current bottlenecks**:

- Uses `deepcopy` for state management
- Iterates over all neighbors repeatedly
- Allocates temporary arrays

**Optimization ideas**:

- Pre-allocate workspace buffers
- Use in-place operations
- Cache neighbor lists
- Reduce `deepcopy` usage

### 2. `get_edges_in_groups` Function

**Location**: `src/assignment.jl`

**Why it matters**: Called during log-likelihood recomputation

**Current bottlenecks**:

- Uses `findall` (allocates)
- Creates new vector each time
- Linear search through nodes

**Optimization ideas**:

- Pre-compute group membership indices
- Use pre-allocated output buffers
- Cache results for frequently accessed group pairs

### 3. Log-likelihood Updates

**Location**: `src/optimization/swap_workspace.jl`, `swap_categorical.jl`

**Why it matters**: Must be computed after each swap

**Current approach**: Recomputes only affected group pairs (good!)

**Optimization ideas**:

- Batch `logpdf` computations
- Use vectorized operations
- Cache intermediate calculations

## 📁 File Structure

```
dev/
  ├── run_benchmarks.jl          # Easy-to-use benchmark runner
  ├── benchmark_optimization.jl   # Standalone benchmarking script
  ├── BENCHMARKING.md            # Detailed documentation
  └── benchmark_results/         # Stored benchmark results
      └── baseline.json          # Reference baseline

test/
  └── test_performance_regression.jl  # Performance tests for CI
```

## 📈 Example Output

```
--- Single Swap Operations (Bernoulli) ---
Benchmarking Bernoulli swap (n=50, k=2)...
  Median: 0.234 ms
Benchmarking Bernoulli swap (n=200, k=3)...
  Median: 1.567 ms

========================================
Performance Comparison vs Baseline
Baseline: 2024-10-15 14:30:00
========================================
✓ FASTER bernoulli_swap_n50_k2: 1.23x (23.0%)
         Current: 0.190 ms | Baseline: 0.234 ms

≈ SIMILAR bernoulli_swap_n200_k3: 1.02x (2.0%)
         Current: 1.537 ms | Baseline: 1.567 ms
```

## 🔧 Advanced Usage

### Run Only Specific Benchmarks

Edit `dev/benchmark_optimization.jl` to comment out benchmarks you don't need.

### Compare Two Specific Benchmark Files

```bash
julia dev/run_benchmarks.jl compare results/v1.json results/v2.json
```

### Profile Your Code

```julia
using Profile

include("dev/test_decorated_paper.jl")

# Profile a specific function
@profile main(500:500:1000, 2)

Profile.print()
# Or for a flamegraph:
using ProfileView
ProfileView.view()
```

### Check Allocations

```julia
using BenchmarkTools

# See allocations for a single operation
@btime apply_swap!($assignment, $swap) samples=1 evals=1
```

## 🎓 Best Practices

1. **Always establish a baseline first** - You need a reference point
2. **Make incremental changes** - Change one thing at a time
3. **Profile before optimizing** - Don't guess where the bottleneck is
4. **Test correctness** - Fast but wrong is useless
5. **Document your changes** - Explain why you made each optimization
6. **Consider maintainability** - Don't sacrifice readability for tiny gains

## 📚 Resources

- **Detailed benchmarking guide**: See `dev/BENCHMARKING.md`
- **Julia Performance Tips**:
  https://docs.julialang.org/en/v1/manual/performance-tips/
- **BenchmarkTools.jl**: https://juliaci.github.io/BenchmarkTools.jl/stable/
- **Profile module**: https://docs.julialang.org/en/v1/stdlib/Profile/

## 🤝 Contributing Performance Improvements

When submitting a PR with performance improvements:

1. Include before/after benchmark results
2. Explain what you optimized and why
3. Ensure all tests still pass
4. Document any trade-offs made
5. Consider adding new benchmarks for your changes

## ❓ Troubleshooting

### "BenchmarkTools not available"

```bash
julia --project=test -e 'using Pkg; Pkg.add("BenchmarkTools")'
```

### High variance in results

- Close other applications
- Run benchmarks multiple times
- Use `--threads=1` flag for consistency

### Benchmark takes too long

- Reduce the `samples` parameter
- Use smaller test networks
- Run individual benchmark categories instead of all at once

## 📞 Getting Help

- Open an issue with benchmark results
- Include your system specs (OS, Julia version, CPU)
- Describe what you're trying to optimize

---

Happy optimizing! 🚀
