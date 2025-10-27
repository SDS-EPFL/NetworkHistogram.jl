module NetworkHistogram
using Accessors
using StatsBase
using StaticArrays
using ProgressMeter
import StatsAPI: loglikelihood, fit, params
import Base: convert, eltype, zero
using Distributions
using LinearAlgebra
using ArgCheck
import Random: randperm, AbstractRNG, rand
import Distributions: logpdf, pdf

using IntervalSets

using Reexport
@reexport using Graphons

import Graphons: _extract_param, convert_to_params

include("utils/include.jl")

@reexport using .FastSymArray

include("distributions/include.jl")
include("EdgeList.jl")
include("assignment.jl")
include("optimization/greedy.jl")
include("distributions/hist_dist.jl")
include("preprocessor/abstractConvertor.jl")
include("preprocessor/categorical.jl")
include("preprocessor/continuous.jl")
include("estimator/abstractEstimator.jl")
include("estimator/SpectralEstimator.jl")
include("api.jl")

export EdgeList, neighbors, nodes, loglikelihood, zero, fit, agg_params, logpdf, pdf,
       GreedyParams, nethist, nethist_discrete_edges, ordered_start_labels, RandomGroupSwap,
       Strict, PreviousBestValue, nethist_binary_edges

function from_adjs_to_decorated end

function heatmap_params end

export from_adjs_to_decorated, heatmap_params

export NethistResult

end
