
# helpers for optimal transport alignment

function plan_to_permutation(plan)
    ordering = argmax(plan, dims = 1) .|> Tuple |> vec
    perm = sort(ordering, by = x -> x[1]) .|> last
    return perm
end

"""
Get the permutation aligning source and target matrices using optimal transport.

This function converts a gromov-wasserstein plan into a permutation by taking the argmax
along the rows.

This function uses [`gromov_wasserstein`](https://pythonot.github.io/gen_modules/ot.gromov.html#ot.gromov.BAPG_gromov_wasserstein)

# See also
- [`align_matrices`](@ref)
- [`ot.gromov.gromov_wasserstein`](@extref)
"""
function get_perm_alignment(
        src::AbstractMatrix{<:Real},
        target::AbstractMatrix{<:Real};
        kwargs...
)
    plan = ot[].gromov.gromov_wasserstein(
        C2 = jl_to_np(src), C1 = jl_to_np(target), kwargs...)
    plan = pyconvert(typeof(target), plan)
    return plan_to_permutation(plan)
end

function get_perm_alignment(
        src::AbstractMatrix{T1},
        target::AbstractMatrix{T2};
        kwargs...
) where {T1 <: AbstractVector, T2 <: AbstractVector}
    C1 = jl_to_np(target)
    C2 = jl_to_np(src)
    dist, log_ = fngw.x.fused_network_gromov_wasserstein2(
        M = jl_to_np(zeros(size(target, 1), size(src, 1))),
        C1 = C1,
        C2 = C2,
        A1 = jl_to_np(ones(size(target, 1), size(target, 1))),
        A2 = jl_to_np(ones(size(src, 1), size(src, 1))),
        p = jl_to_np(fill(1.0 / size(target, 1), size(target, 1))),
        q = jl_to_np(fill(1.0 / size(src, 1), size(src, 1))),
        alpha = 1.0,
        beta = 0.0,
        log = true,
        kwargs...
    )
    plan = pyconvert(Matrix{Float64}, log_["T"])
    return plan_to_permutation(plan)
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
