"""
    SBMEstimator

Abstract base type for all Stochastic Block Model (SBM) estimators.

All concrete estimator types should implement:
- `estimate(estimator, data, initial_labels; progress=true)`: Main estimation function
- `score(estimator)`: Return current objective value (if applicable)
"""
abstract type SBMEstimator end

abstract type Result end

# struct NethistResult{L, M} <: Result
#     labels::L
#     model::M
# end

include("GreedySuffStats.jl")
