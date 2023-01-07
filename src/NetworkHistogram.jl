module NetworkHistogram

using ValueHistories

include("assignment.jl")
include("optimize.jl")
include("histogram.jl")

include("update/accept_reject.jl")
include("update/select_swap.jl")
include("update/proposal.jl")

end
