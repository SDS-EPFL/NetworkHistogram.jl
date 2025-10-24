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
using .FastSymArray

include("distributions/include.jl")
include("EdgeList.jl")
include("assignment.jl")
include("optimization/greedy.jl")
include("estimator/abstractEstimator.jl")
include("estimator/SpectralEstimator.jl")
include("api.jl")

export EdgeList, neighbors, nodes, loglikelihood, zero, fit, agg_params, logpdf

function from_adjs_to_decorated end

function heatmap_params end

export from_adjs_to_decorated, heatmap_params
end
