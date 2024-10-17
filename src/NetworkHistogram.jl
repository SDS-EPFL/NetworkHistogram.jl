module NetworkHistogram

using LinearAlgebra, SparseArrays, DataStructures
using Distributions, DensityInterface
using Graphs
using PermutationSymmetricTensors
using ProgressMeter: Progress, next!, finish!
import StatsBase, Random
using DensityInterface: logdensityof

import Distributions.fit

include("assignments/Assignments.jl")
include("sbm.jl")
include("observations.jl")
include("fit.jl")
include("optimisation/include.jl")

# more specialised and faster assignment types and methods
include("assignments/include.jl")

@warn "User interface is not yet implemented"

end
