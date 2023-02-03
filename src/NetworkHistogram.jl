module NetworkHistogram

using ValueHistories

include("group_numbering.jl")
include("assignment.jl")
include("optimize.jl")
include("histogram.jl")
include("proposal/select_swap.jl")
include("proposal/proposal.jl")

end
