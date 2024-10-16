module NetworkHistogram

using LinearAlgebra, SparseArrays, DataStructures
using Distributions, DensityInterface
using Graphs
using PermutationSymmetricTensors
using ProgressMeter: Progress,next!,finish!
import StatsBase, Random

include("group_numbering.jl")
include("sbm.jl")
include("observations.jl")
include("assignments/Assignments.jl")
include("fit.jl")
include("optimisation/include.jl")

# more specialised and faster assignment types and methods
include("assignments/include.jl")

end
