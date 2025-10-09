module NetworkHistogram
using StatsBase
using StaticArrays
using ProgressMeter
import StatsAPI: loglikelihood, fit, params
import Base: convert, eltype, zero
using Distributions

include("utils/include.jl")
using .FastSymArray

include("distributions/include.jl")
include("EdgeList.jl")
include("assignment.jl")
include("block_model.jl")
include("optimization/greedy.jl")
include("api.jl")

export EdgeList, neighbors, nodes, loglikelihood, zero, fit, agg_params, logpdf

function from_adjs_to_decorated end

function heatmap_params end

export from_adjs_to_decorated, heatmap_params
end
