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

include("distributions/hist_dist.jl")
include("preprocessor/abstractConvertor.jl")
include("estimator/abstractEstimator.jl")
include("api.jl")

export GreedyParams, nethist, nethist_discrete_edges, ordered_start_labels, RandomGroupSwap,
       Strict, PreviousBestValue, nethist_binary_edges

end
