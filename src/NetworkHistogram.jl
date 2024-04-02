module NetworkHistogram

using ValueHistories, StatsBase, Random, LinearAlgebra, Kronecker
using Clustering: kmeans

using Arpack: eigs
using ArnoldiMethod: partialschur, partialeigen, SR, LR

export graphhist, PreviousBestValue, Strict, RandomNodeSwap
export OrderedStart, RandomStart, EigenStart, DistStart

include("group_numbering.jl")
include("assignment.jl")
include("history.jl")

include("config_rules/starting_assignment_rule.jl")
include("config_rules/swap_rule.jl")
include("config_rules/accept_rule.jl")
include("config_rules/stop_rule.jl")
include("config_rules/bandwidth_selection_rule.jl")

include("optimize.jl")
include("histogram.jl")
include("proposal.jl")

include("utils/adjacency_utils.jl")
include("utils/utils.jl")

include("data/gt.jl")
include("data/datasets.jl")
end
