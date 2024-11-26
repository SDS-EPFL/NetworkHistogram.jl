module NetworkHistogram

using LinearAlgebra, SparseArrays
using Distributions, DensityInterface
using Graphs, SimpleWeightedGraphs
using PermutationSymmetricTensors
using ProgressMeter: Progress, next!, finish!
import StatsBase, Random
using DensityInterface: logdensityof
using StaticArrays: MVector, MMatrix
using LogExpFunctions: xlogx, xlogy
using ArnoldiMethod: LM, SR, LR, partialschur, partialeigen
import Arpack
import Metis
import IterativeSolvers
import Clustering
import StatsAPI: loglikelihood, fit
using CategoricalArrays, CategoricalDistributions
using Discretizers: LinearDiscretizer, binedges, DiscretizeUniformWidth, encode
using Combinatorics: permutations

include("assignments/Assignments.jl")
include("sbm.jl")
include("observations.jl")
include("optimisation/include.jl")

# more specialised and faster assignment types and methods
include("assignments/include.jl")

@warn "User interface is not yet implemented"

export loglikelihood, fit

# export options for optimisation
export estimate_graphon
# starting assignment rules
export InitRule
export OrderedStart, RandomStart, SpectralStart, MetisStart, FromAssignment
# accept rules
export AcceptRule
export Strict
# stopping rules
export PreviousBestValue
# bandwidth selection rules
export OracleK, EstimatedEigenvalues, EstimatedDegrees, select_number_node_per_block
# random local search rules
export RandomNodeSwap, RandomGroupSwap

# export useful function for manipulating assignments
export Assignment, number_groups, number_nodes
export get_ordered_adjacency_matrix, get_vertex_in_group, get_group_of_vertex
export BernoulliData, CategoricalData
export Observations, discretise

end
