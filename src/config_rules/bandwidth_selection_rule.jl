
function select_bandwidth(A::Array{T, 2}; type = "degs", alpha = 1)::Int where {T}
    h = oracle_bandwidth(A; type=type, alpha=alpha)
    return max(2, min(size(A)[1], round(Int, h)))
end

function select_bandwidth(A::Array{T, 3}; type = "degs", alpha = 1)::Int where {T}
    hs = [select_bandwidth(A[:, :, i]; type, alpha) for i in 1:size(A, 3)]
    #@warn "Naive bandwidth selection for multilayer graph histogram: using minimum over layers"
    h = max(2, min(size(A, 1), round(Int, minimum(hs))))
    return h
end

"""
    oracle_bandwidth(A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))

Oracle bandwidth selection for graph histogram, using

```math
\\widehat{h^*}=\\left(2\\left(\\left(d^T d\\right)^{+}\\right)^2 d^T A d \\cdot \\hat{m} \\hat{b}\\right)^{-\\frac{1}{2}} \\hat{\\rho}_n^{\\frac{1}{4}},
```

where ``d`` is the vector of degree sorted in increasing order,``\\hat{\\rho}_n`` is the empirical edge density, and  ``m``, ``b`` are the slope and intercept fitted on ``d[n/2-c\\sqrt{n}:n/2+c\\sqrt{n}]`` for some ``c``.
"""
function oracle_bandwidth(A; type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8), diff_degrees = 5)
    if type ∉ ["eigs", "degs"]
        error("Invalid input type $(type)")
    end
    A = drop_disconnected_components(A)

    n = size(A, 1)
    bound =  round(Int, c * sqrt(n))
    rho = sum(A) / (n * (n - 1))
    rhoHat_inv = inv(rho)

    # Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
    if type == "eigs"
        eig_res = eigs(A, nev = 1, which = :LM)
        u = eig_res[2]
        mult = eig_res[1][1]
    elseif type == "degs"
        u = sum(A, dims = 2)
        mult = (u' * A * u) / (sum(u .^ 2))^2
    else
        error("Invalid input type $(type)")
    end

    # Calculation bandwidth
    u = sort(u, dims = 1)
    uMid = u[max(1, n ÷ 2 - bound):min(n, n ÷ 2 + bound)]

    tries = 0
    increment = 1
    while length(unique(uMid)) < diff_degrees && tries < 10 && bound < n÷2
        bound += increment
        uMid = u[max(1, n ÷ 2 - bound):min(n, n ÷ 2 + bound)]
        tries += 1
    end
    js = collect(max(1, n ÷ 2 - bound):min(n, n ÷ 2 + bound)) .- n÷2
    coeff = hcat(ones(length(uMid)), js) \ uMid
    h = (2^(alpha) * alpha^(1/2) * mult * coeff[1] * coeff[2] )^(-1 /(alpha+1))*rho^(1/(2*(alpha+1)))

    return h[1]
end
