abstract type SuffStats end

function add_sample end
function remove_sample end
function make_k_block end

# score will be minimized
function score end
function to_params end

# some suffstat may need the edge index (i,j) to update properly
add_sample(suffstats::SuffStats, sample, i, j) = add_sample(suffstats, sample)
remove_sample(suffstats::SuffStats, sample, i, j) = remove_sample(suffstats, sample)

include("categorical.jl")
include("bernoulli.jl")
include("generic.jl")
