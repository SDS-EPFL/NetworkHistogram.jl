module PythonOptimalTransport
using PythonCall
using NetworkHistogram

import NetworkHistogram: align_matrices, get_perm_alignment

const ot = Ref{Py}()
const fngw = Ref{Py}()

function __init__()
    ot[] = pyimport("ot")
    pyimport("sys").path.append(@__DIR__)
    fngw[] = pyimport("fngw")
end

jl_to_np(mat) = Py(mat).to_numpy()

"""
Get the permutation aligning source and target matrices using optimal transport.

This function converts a gromov-wasserstein plan into a permutation by taking the argmax
along the rows.

This function uses [`gromov_wasserstein`](https://pythonot.github.io/gen_modules/ot.gromov.html#ot.gromov.BAPG_gromov_wasserstein)

# See also
- [`align_matrices`](@ref)
- [`ot.gromov.gromov_wasserstein`](@extref)
"""
function get_perm_alignment(src, target; kwargs...)
    plan = ot[].gromov.gromov_wasserstein(
        C2 = jl_to_np(src), C1 = jl_to_np(target), kwargs...)
    plan = pyconvert(Matrix{Float64}, plan)
    ordering = argmax(plan, dims = 1) .|> Tuple |> vec
    perm = sort(ordering, by = x -> x[1]) .|> last
    return perm
end

"""
Align the source and target matrices using optimal transport.

# See also
- [`get_perm_alignment`](@ref).
"""
function align_matrices(src, target)
    perm = get_perm_alignment(src, target)
    return src[perm, perm]
end

end
