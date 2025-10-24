module NetworkHistogram
using StatsBase
using StaticArrays
using ProgressMeter
import StatsAPI: loglikelihood, fit, params
import Base: convert, eltype, zero
using Distributions
using LinearAlgebra
using ArgCheck
using Random: randperm

using Reexport
@reexport using Graphons

include("utils/include.jl")

@reexport using .FastSymArray

include("distributions/include.jl")
include("EdgeList.jl")
include("assignment.jl")
include("optimization/greedy.jl")
include("preprocessor/categorical.jl")
include("preprocessor/continuous.jl")
include("estimator/abstractEstimator.jl")
include("estimator/SpectralEstimator.jl")
include("api.jl")

export EdgeList, neighbors, nodes, loglikelihood, zero, fit, agg_params, logpdf,
       GreedyParams, nethist, nethist_discrete_edges, ordered_start_labels, RandomGroupSwap,
       Strict, PreviousBestValue, nethist_binary_edges

function from_adjs_to_decorated end

function heatmap_params end

export from_adjs_to_decorated, heatmap_params

export NethistResult

end
