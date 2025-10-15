# Performance Regression Test Suite - Summary

## What Was Created

A comprehensive performance benchmarking and profiling suite for
NetworkHistogram optimization, consisting of:

### 1. **Test Files**

- `test/test_performance_regression.jl` - Performance regression tests that run
  with `Pkg.test()`
  - Tests for Bernoulli and Categorical networks
  - Multiple network sizes (50, 200, 500 nodes)
  - Single swap operations and full optimization workflows
  - Component-level benchmarks

### 2. **Standalone Benchmarking**

- `benchmark_optimization.jl` - Comprehensive standalone benchmark suite
  - Saves results to JSON with timestamps
  - Automatic comparison with baseline
  - Detailed performance metrics (median, mean, std, min, max)

### 3. **Easy-to-Use Runner**

- `dev/run_benchmarks.jl` - User-friendly command-line interface
  - Simple commands: `baseline`, `current`, `compare`, `clean`
  - Handles dependencies automatically
  - Interactive confirmations for destructive operations

### 4. **Profiling Tools**

- `dev/profile_optimization.jl` - Profiling helper
  - Profile swap operations, full optimization, or components
  - Integrated flamegraph support
  - Configurable network sizes and iteration counts

### 5. **Documentation**

- `PERFORMANCE.md` - Main performance guide
- `benchmark/BENCHMARKING.md` - Detailed benchmarking documentation
- This summary document

## How to Use

### Quick Start (5 minutes)

```bash
# 1. Create baseline
julia dev/run_benchmarks.jl baseline

# 2. Make your optimizations in src/optimization/

# 3. Test performance
julia dev/run_benchmarks.jl current

# 4. Verify correctness
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Example Output

```
--- Single Swap Operations (Bernoulli) ---
Benchmarking Bernoulli swap (n=50, k=2)...
  Median: 0.234 ms

========================================
Performance Comparison vs Baseline
========================================
✓ FASTER bernoulli_swap_n50_k2: 1.23x (23.0%)
         Current: 0.190 ms | Baseline: 0.234 ms
```

## Key Insights from Code Analysis

Based on analysis of `test_decorated_paper.jl` and the source code:

### Primary Bottlenecks

1. **`apply_swap!`** (called millions of times)

   - Location: `src/optimization/swap_workspace.jl`, `swap_categorical.jl`
   - Issues: Uses `deepcopy`, allocates during iteration
   - Impact: 🔴 CRITICAL - dominates runtime

2. **`get_edges_in_groups`** (called during LL updates)

   - Location: `src/assignment.jl`
   - Issues: Uses `findall`, allocates new vectors
   - Impact: 🟡 MODERATE - called less frequently

3. **Edge iteration** (used throughout)
   - Location: `src/EdgeList.jl`
   - Issues: Iterator overhead
   - Impact: 🟢 LOW - but cumulative

### Workflow from test_decorated_paper.jl

The typical optimization workflow:

1. Create SBM (Stochastic Block Model)
2. Sample network from SBM
3. Initialize node labels
4. Run greedy optimization with `nethist()`
   - Iteratively swap nodes between groups
   - Accept swaps that improve log-likelihood
5. Measure convergence via log-likelihood

## Benchmarked Scenarios

### Network Sizes

- **Small**: n=50 nodes, k=2 groups (quick iteration)
- **Medium**: n=200 nodes, k=3 groups (realistic size)
- **Large**: n=500 nodes, k=5 groups (stress test)

### Network Types

- **Bernoulli**: Binary edges (0/1) - simpler, faster
- **Categorical**: Multi-valued edges (m categories) - more complex

### Benchmark Types

- **Single swap**: Apply + revert one node swap
- **Full optimization**: Complete optimization run (1k iterations)
- **Components**: Individual function benchmarks

## Files and Their Purpose

```
NetworkHistogram/
├── PERFORMANCE.md                    # Main guide (START HERE)
├── test/
│   ├── test_performance_regression.jl  # CI-friendly tests
│   └── Project.toml                   # Added BenchmarkTools dependency
├── dev/
│   ├── run_benchmarks.jl              # 👈 Easy CLI (USE THIS)
│   ├── benchmark_optimization.jl      # Core benchmarking logic
│   ├── profile_optimization.jl        # Profiling helper
│   ├── BENCHMARKING.md               # Detailed docs
│   ├── benchmark_results/            # Stored results (auto-created)
│   │   └── baseline.json             # Your reference baseline
│   └── test_decorated_paper.jl       # Original workflow example
└── src/optimization/                 # 🎯 Optimize these files
    ├── greedy.jl
    ├── swap_workspace.jl
    ├── swap_categorical.jl
    └── config_rules/
```

## Common Workflows

### A. Making Performance Improvements

```bash
# Step 1: Baseline
julia dev/run_benchmarks.jl baseline

# Step 2: Profile to find bottlenecks
julia dev/profile_optimization.jl swap

# Step 3: Make changes to src/optimization/

# Step 4: Benchmark
julia dev/run_benchmarks.jl current

# Step 5: Test correctness
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 6: Repeat steps 2-5 until satisfied
```

### B. Comparing Two Versions

```bash
# Benchmark version A
git checkout feature-A
julia dev/run_benchmarks.jl results_A.json

# Benchmark version B
git checkout feature-B
julia dev/run_benchmarks.jl results_B.json

# Compare
julia dev/run_benchmarks.jl compare results_A.json results_B.json
```

### C. Debugging Performance Regression

```bash
# Find when regression occurred
git bisect start
git bisect bad HEAD
git bisect good v1.0.0

# For each commit
julia dev/run_benchmarks.jl
# Mark good/bad based on results
git bisect good  # or bad
```

## Optimization Strategies

### 1. Profile First

Don't guess - use `profile_optimization.jl` to see what's actually slow.

### 2. Reduce Allocations

The biggest wins usually come from eliminating allocations in hot paths.

**Check allocations**:

```julia
using BenchmarkTools
@btime apply_swap!($assignment, $swap) samples=1 evals=1
#           ^^^^^ This shows allocations
```

**Common fixes**:

- Pre-allocate buffers
- Use `@inbounds` (after bounds checking once)
- Avoid `deepcopy` when possible
- Use views instead of copies

### 3. Type Stability

Julia is fast when types are known at compile time.

**Check type stability**:

```julia
using Cthulhu
@descend apply_swap!(assignment, swap)
# Look for red (runtime dispatch)
```

### 4. SIMD/Vectorization

For bulk operations on arrays, help the compiler vectorize.

### 5. Cache-Friendly Access

Access memory in order when possible (column-major for Julia).

## Expected Performance Gains

Based on typical optimization opportunities in similar codebases:

- **Low-hanging fruit** (reduce allocations): 20-50% speedup
- **Algorithm improvements** (better data structures): 2-10x speedup
- **SIMD/vectorization**: 2-4x speedup (for vectorizable operations)
- **Type stability fixes**: 2-5x speedup (if unstable)

The swap operation is called O(iterations × n) times, so even small
improvements compound significantly.

## Integration with CI/CD

Add to `.github/workflows/benchmark.yml`:

```yaml
name: Benchmark
on: [pull_request]

jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2

      - uses: julia-actions/setup-julia@v1

      - name: Benchmark PR
        run: julia benchmark_optimization.jl pr_results.json

      - name: Benchmark main
        run: |
          git fetch origin main
          git checkout origin/main
          julia benchmark_optimization.jl main_results.json

      - name: Compare
        run:
          julia dev/run_benchmarks.jl compare pr_results.json main_results.json
```

## Troubleshooting

### "BenchmarkTools not found"

```bash
julia --project=test -e 'using Pkg; Pkg.add("BenchmarkTools")'
```

### Results vary too much

- Close other applications
- Disable CPU frequency scaling
- Run with `--threads=1`
- Increase sample count

### Benchmark takes too long

- Reduce `samples` parameter
- Use smaller networks
- Run specific benchmarks only

## Next Steps

1. **Establish your baseline**: `julia dev/run_benchmarks.jl baseline`
2. **Read the detailed docs**: See `benchmark/BENCHMARKING.md`
3. **Profile the code**: `julia dev/profile_optimization.jl swap`
4. **Start optimizing**: Focus on `apply_swap!` first
5. **Measure improvements**: `julia dev/run_benchmarks.jl current`
6. **Share results**: Open PR with before/after benchmarks

## Questions?

- Check `PERFORMANCE.md` for main guide
- Check `benchmark/BENCHMARKING.md` for detailed docs
- Run `julia dev/run_benchmarks.jl help`
- Run `julia dev/profile_optimization.jl help`

---

Happy optimizing! The suite is designed to make performance work systematic and
data-driven. 🚀
