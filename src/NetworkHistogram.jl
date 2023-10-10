module NetworkHistogram

using ValueHistories, StatsBase, Random, LinearAlgebra, Kronecker

using Arpack: eigs

export graphhist, PreviousBestValue, Strict, RandomNodeSwap
export OrderedStart, RandomStart, EigenStart

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

include("utils.jl")

include("data/gt.jl")
include("data/datasets.jl")
end
