@warn "Deprecated bandwidth selection needs to be updated"

function select_bandwidth(A::Array{T, 2}; type = "degs", alpha = 1, c = 1)::Int where {T}
    h = oracle_bandwidth(A, type, alpha, c)
    return max(2, min(size(A)[1], round(Int, h)))
end

"""
    oracle_bandwidth(A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))

Oracle bandwidth selection for graph histogram, using

```math
\\widehat{h^*}=\\left(2\\left(\\left(d^T d\\right)^{+}\\right)^2 d^T A d \\cdot \\hat{m} \\hat{b}\\right)^{-\\frac{1}{2}} \\hat{\\rho}_n^{\\frac{1}{4}},
```

where ``d`` is the vector of degree sorted in increasing order,``\\hat{\\rho}_n`` is the empirical edge density, and  ``m``, ``b`` are the slope and intercept fitted on ``d[n/2-c\\sqrt{n}:n/2+c\\sqrt{n}]`` for some ``c``.
"""
function oracle_bandwidth(A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))
    if type ∉ ["eigs", "degs"]
        error("Invalid input type $(type)")
    end

    if alpha != 1
        error("Currently only supports alpha = 1")
    end

    n = size(A, 1)
    midPt = collect(max(1, round(Int, (n ÷ 2 - c * sqrt(n)))):round(Int,
        (n ÷ 2 + c * sqrt(n))))
    rhoHat_inv = inv(sum(A) / (n * (n - 1)))

    # Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
    if type == "eigs"
        eig_res = eigs(A, nev = 1, which = :LM)
        u = eig_res.vectors
        mult = eig_res.values[1]
    elseif type == "degs"
        u = sum(A, dims = 2)
        mult = (u' * A * u) / (sum(u .^ 2))^2
    else
        error("Invalid input type $(type)")
    end

    # Calculation bandwidth
    u = sort(u, dims = 1)
    uMid = u[midPt]
    lmfit_coef = hcat(ones(length(uMid)), 1:length(uMid)) \ uMid

    h = (2^(alpha + 1) * alpha * mult^2 *
         (lmfit_coef[2] * length(uMid) / 2 + lmfit_coef[1])^2 * lmfit_coef[2]^2 *
         rhoHat_inv)^(-1 / (2 * (alpha + 1)))
    #estMSqrd = 2*mult^2*(lmfit_coef[2]*length(uMid)/2+lmfit_coef[1])^2*lmfit_coef[2]^2*rhoHat_inv^2*(n+1)^2
    return h[1]
end
