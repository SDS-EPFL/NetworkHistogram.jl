module NetworkHistogram
using StatsBase
using StaticArrays

include("utils/include.jl")
using .FastSymArray

include("distributions_type.jl")
include("block_model.jl")
include("EdgeList.jl")
include("assignment.jl")
include("optimization/greedy.jl")


export EdgeList, neighbors, nodes

end
