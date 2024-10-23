module NetworkHistogram

using LinearAlgebra, SparseArrays, DataStructures
using Distributions, DensityInterface
using Graphs, SimpleWeightedGraphs
using PermutationSymmetricTensors
using ProgressMeter: Progress, next!, finish!
import StatsBase, Random
using DensityInterface: logdensityof
using StaticArrays: MVector, MMatrix
using LogExpFunctions: xlogx, xlogy
using LoopVectorization: @turbo
using ArnoldiMethod: LM, SR, LR, partialschur, partialeigen
import Arpack
import Metis
import IterativeSolvers
import Clustering
import StatsAPI: loglikelihood, fit
using CategoricalArrays, CategoricalDistributions
export loglikelihood, fit

include("assignments/Assignments.jl")
include("sbm.jl")
include("observations.jl")
include("optimisation/include.jl")

# more specialised and faster assignment types and methods
include("assignments/include.jl")

@warn "User interface is not yet implemented"

end
