module BootstrapExt

using NetworkHistogram
import NetworkHistogram: test_extension_boot

using Bootstrap

function test_extension_boot()
    return "Bootstrap extension works!"
end

end
