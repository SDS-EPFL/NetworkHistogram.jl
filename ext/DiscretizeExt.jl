module DiscretizeExt

using NetworkHistogram
import NetworkHistogram: test_extension_disc, get_ref_dist, Dist, ZeroInflated
import Distributions: ContinuousUnivariateDistribution
using DiscretizeDistributions

# in_interval.(x, support(discretized_dist))
function test_extension_disc()
    return "Discretize extension works!"
end

function get_ref_dist(dist::D, ::Val{true}) where {D <: ContinuousUnivariateDistribution}
    return Dist(ZeroInflated(dist))
end
function get_ref_dist(dist::D, ::Val{false}) where {D <: ContinuousUnivariateDistribution}
    return Dist(dist)
end

end
