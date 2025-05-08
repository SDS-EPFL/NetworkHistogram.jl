module NetworkHistogram
using StatsBase
using StaticArrays
using ProgressMeter
import StatsAPI: loglikelihood
import Base: convert, eltype

include("utils/include.jl")
using .FastSymArray

include("distributions_type.jl")
include("block_model.jl")
include("EdgeList.jl")
include("assignment.jl")
include("optimization/greedy.jl")
include("api.jl")

export EdgeList, neighbors, nodes, loglikelihood

end
