# Performance Optimization Quick Reference

## 🚀 Quick Commands

```bash
# Establish baseline
julia dev/run_benchmarks.jl baseline

# Benchmark current code
julia dev/run_benchmarks.jl current

# Profile to find bottlenecks
julia dev/profile_optimization.jl swap

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Visualize results
julia dev/visualize_benchmarks.jl --all
```

## 📊 Understanding Output

```
✓ FASTER   = >5% improvement
✗ SLOWER   = >5% regression
≈ SIMILAR  = Within ±5%
```

## 🎯 Priority Hotspots

### 1. `apply_swap!` 🔴 CRITICAL
- **File**: `src/optimization/swap_workspace.jl`, `swap_categorical.jl`
- **Why**: Called ~1M times per run
- **Fix**: Reduce allocations, avoid `deepcopy`

### 2. `get_edges_in_groups` 🟡 MODERATE
- **File**: `src/assignment.jl`
- **Why**: Called during LL updates
- **Fix**: Pre-allocate, cache group membership

### 3. Edge iteration 🟢 LOW
- **File**: `src/EdgeList.jl`
- **Why**: Used everywhere
- **Fix**: Ensure type stability

## 🛠️ Common Optimizations

### Check Allocations
```julia
using BenchmarkTools
@btime my_function($args) samples=1 evals=1
# Look for allocations in output
```

### Profile Code
```julia
using Profile
@profile my_function(args)
Profile.print(maxdepth=15)
```

### Type Stability
```julia
using Cthulhu
@descend my_function(args)
# Red = type unstable (BAD)
```

## 📁 Key Files

```
├── PERFORMANCE.md               # Main guide
├── dev/
│   ├── run_benchmarks.jl       # 👈 USE THIS
│   ├── profile_optimization.jl  # For profiling
│   ├── visualize_benchmarks.jl  # View results
│   └── BENCHMARKING.md         # Details
└── src/optimization/           # 🎯 Optimize here
    ├── swap_workspace.jl
    └── swap_categorical.jl
```

## 📈 Expected Gains

- Reduce allocations: **20-50%** speedup
- Better data structures: **2-10x** speedup
- SIMD/vectorization: **2-4x** speedup
- Fix type instability: **2-5x** speedup

## 🔄 Workflow

1. **Baseline** → 2. **Profile** → 3. **Optimize** → 4. **Benchmark** → 5. **Test** → Repeat

## 💡 Tips

- Focus on **hot paths** (profile first!)
- Measure **before and after** every change
- Keep changes **small and focused**
- Always **test correctness**
- Document **what and why**

## 🆘 Troubleshooting

### "BenchmarkTools not found"
```bash
julia --project=test -e 'using Pkg; Pkg.add("BenchmarkTools")'
```

### Results vary
- Close other apps
- Use `--threads=1`
- Increase samples

### Too slow
- Reduce samples
- Use smaller networks
- Run specific benchmarks

## 📚 Learn More

- `PERFORMANCE.md` - Full guide
- `dev/BENCHMARKING.md` - Detailed docs
- `julia dev/run_benchmarks.jl help` - CLI help

---

**Remember**: Profile → Optimize → Benchmark → Test → Repeat 🔁
