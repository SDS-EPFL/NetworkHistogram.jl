module DiscretizeExt

using NetworkHistogram
import NetworkHistogram: get_ref_dist, Dist, ZeroInflated
import Distributions: ContinuousUnivariateDistribution
using DiscretizeDistributions

function get_ref_dist(dist::D, ::Val{true}) where {D <: ContinuousUnivariateDistribution}
    return Dist(ZeroInflated(dist))
end
function get_ref_dist(dist::D, ::Val{false}) where {D <: ContinuousUnivariateDistribution}
    return Dist(dist)
end

end
