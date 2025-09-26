module DistExt

using NetworkHistogram
import NetworkHistogram: test_extension_dist

using Distributions

function test_extension_dist()
    return "Distribution extension works!"
end

end
