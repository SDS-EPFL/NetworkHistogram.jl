module NetworkHistogram

using LinearAlgebra, SparseArrays, DataStructures
using Distributions, DensityInterface
using Graphs
using PermutationSymmetricTensors

import StatsBase, Random

include("group_numbering.jl")
include("sbm.jl")
include("observations.jl")
include("assignments/Assignments.jl")
include("optimisation/swap.jl")
include("optimisation/opti.jl")
include("assignments/data_structures.jl")

end
