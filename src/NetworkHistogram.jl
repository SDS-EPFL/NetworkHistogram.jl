module NetworkHistogram

using LinearAlgebra, SparseArrays
using Distributions, DensityInterface
using Graphs, SimpleWeightedGraphs
using PermutationSymmetricTensors
using ProgressMeter: Progress, next!, finish!, ProgressUnknown
import StatsBase, Random
using DensityInterface: logdensityof
using LogExpFunctions: xlogx, xlogy
using ArnoldiMethod: LM, SR, LR, partialschur, partialeigen
using KrylovKit: eigsolve
import Metis
import IterativeSolvers
import Clustering
import StatsAPI: loglikelihood, fit
using CategoricalArrays, CategoricalDistributions
using Combinatorics: permutations
using StaticArrays
using Bootstrap: BootstrapSampling, ParametricBootstrapSample, tx, nrun,
                 zeros_tuple
import Bootstrap: bootstrap
import Base.maximum, Base.minimum
import Random: rand
import Base.convert
import Distributions: pdf, logpdf, ncategories, cdf, rand

include("distributions/include.jl")
include("assignments/Assignments.jl")
include("sbm.jl")
include("observations.jl")
include("optimisation/include.jl")

# more specialised and faster assignment types and methods
include("assignments/include.jl")

include("api.jl")
include("bootstrap.jl")

export nethist, nethist_discretised
export loglikelihood, fit, cdf, pdf

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
export OracleK, EstimatedEigenvalues, EstimatedDegrees,
       select_number_node_per_block
# random local search rules
export RandomNodeSwap, RandomGroupSwap

# export useful function for manipulating assignments
export Assignment, number_groups, number_nodes
export get_ordered_adjacency_matrix, get_vertex_in_group, get_group_of_vertex
export BernoulliData, CategoricalData
export Observations, discretise
export DiscretizedDistribution

export Observations, estimate_graphon, nethist, nethist_discretised

export bootstrap

end
