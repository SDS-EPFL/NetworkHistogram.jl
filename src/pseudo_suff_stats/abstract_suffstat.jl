abstract type SuffStats end

function add_sample end
function remove_sample end

add_sample(suffstats::SuffStats, sample, i, j) = add_sample(suffstats, sample)
remove_sample(suffstats::SuffStats, sample, i, j) = remove_sample(suffstats, sample)
function make_k_block end
function score end

include("categorical.jl")
include("bernoulli.jl")
include("generic.jl")
