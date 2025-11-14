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
import Random: randperm, AbstractRNG, rand, shuffle
import Distributions: logpdf, pdf
import LogExpFunctions: xlogx
using IntervalSets
using Hungarian
using Reexport
@reexport using Graphons

import Graphons: _extract_param, convert_to_params, node_labels_to_latents

include("SymArray.jl")
@reexport using .FastSymArray

include("distributions/hist_dist.jl")
include("preprocessor/abstractConvertor.jl")
include("config_rules/include.jl")
include("pseudo_suff_stats/abstract_suffstat.jl")
include("GreedySuffStats.jl")
include("utils/utils_node_labels.jl")
include("api.jl")

export GreedyParams, nethist, nethist_discrete_edges, ordered_start_labels, RandomGroupSwap,
       Strict, PreviousBestValue, nethist_binary_edges

end
