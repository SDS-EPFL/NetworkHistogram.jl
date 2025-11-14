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

function get_perm_alignment(src, target)
    plan = ot[].gromov.gromov_wasserstein(
        C2 = jl_to_np(src), C1 = jl_to_np(target))
    plan = pyconvert(Matrix{Float64}, plan)
    ordering = argmax(plan, dims = 1) .|> Tuple |> vec
    perm = sort(ordering, by = x -> x[1]) .|> last
    return perm
end

function align_matrices(src, target)
    perm = get_perm_alignment(src, target)
    return src[perm, perm]
end

end
