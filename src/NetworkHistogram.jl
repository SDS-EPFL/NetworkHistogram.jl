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

function test_extension_dist end
function test_extension_boot end
function test_extension_disc end

export test_extension_dist, test_extension_boot, test_extension_disc

end
