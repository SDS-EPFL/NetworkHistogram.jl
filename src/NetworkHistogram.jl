module NetworkHistogram

using LinearAlgebra, SparseArrays
using Distributions,DensityInterface
using Graphs
using PermutationSymmetricTensors

import StatsBase, Random

include("group_numbering.jl")
include("sbm.jl")
include("assignments/Assignments.jl")

end
