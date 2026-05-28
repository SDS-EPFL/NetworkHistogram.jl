module PythonOptimalTransport
using PythonCall
using NetworkHistogram

import NetworkHistogram: align_matrices, get_perm_alignment

const ot = Ref{Py}()
const fngw = Ref{Py}()

function __init__()
    ot[] = pyimport("ot")
    pyimport("sys").path.append(@__DIR__)
    # TODO: find why I need to use fngw.x to access the functions later...
    fngw[] = pyimport("fngw")
end

## helpers to convert Julia arrays to numpy arrays
jl_to_np(mat::AbstractArray{<:Real}) = Py(mat).to_numpy()
function jl_to_np(mat::AbstractMatrix{<:AbstractVector})
    Py(permutedims(stack(mat), (3, 2, 1))).to_numpy()
end

include("alignment.jl")

# look at https://pythonot.github.io/auto_examples/backends/plot_optim_gromov_pytorch.html#sphx-glr-auto-examples-backends-plot-optim-gromov-pytorch-py
# to implement semi-relaxed gromov-wasserstein ?
end
