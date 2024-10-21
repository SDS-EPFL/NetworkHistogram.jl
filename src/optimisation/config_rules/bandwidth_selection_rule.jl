abstract type KSelectionRule end
struct OracleK <: KSelectionRule
    K::Int
end
struct OracleM <: KSelectionRule
    M::Int
    α::Float64
end

function OracleM(M::Int)
    return OracleM(M, 1.0)
end

abstract type EstimatedM <: KSelectionRule end
struct EstimatedEigenvalues <: EstimatedM end
struct EstimatedDegrees <: EstimatedM end

"""
    select_number_node_per_block(g::Observations, rule::KSelectionRule)

How to select the number of blocks `K` for the SBM model.

# Implemented rules
- `OracleK(K::Int)`: Use the oracle number of blocks `K`.
- `OracleM(M::Int)`: Give the Holder constant `M` of the graphon, use the results from
    [Olhede and Wolfe (2014)](https://www.pnas.org/doi/epdf/10.1073/pnas.1400374111) to
    estimate the number of blocks `K`.
- `EstimatedEigenvalues()`: Use the estimated eigenvalues of the adjacency matrix to
    estimate the Holder constant and then use `OracleM` to estimate the number of blocks `K`.
- `EstimatedDegrees()`: Use the estimated degrees of the adjacency matrix to estimate the
    Holder constant and then use `OracleM` to estimate the number of blocks `K`.


!!! info
    - The number of blocks `K` should be at most `n/2` where `n` is the number of nodes in
        the graph.
    - The estimated Holder constant `M` comes from equation (11) in Olhede and Wolfe (2014).
"""
select_number_node_per_block

function select_number_node_per_block(g::Observations, rule::OracleK)
    if rule.K > number_nodes(g)÷2
        error("The number of blocks $K is too large for the number of nodes \
        $(number_nodes(g)), it should be at most $(number_nodes(g)÷2)")
    end
    return rule.K
end

function select_number_node_per_block(g::Observations, rule::OracleM)
    rho = density(g.graph)
    n = number_nodes(g)
    k = max(2, round(Int,2*rule.M^2*rho)^(-1/4)*sqrt(n))
    return select_number_node_per_block(g, OracleK(k))
end

function select_number_node_per_block(g::Observations, rule::EstimatedM)
    n = number_nodes(g)
    c = min(4, sqrt(n) / 8)
    number_points_from_mid = round(Int, c * sqrt(n))
    mid_points = collect(max(1, n÷2-number_points_from_mid):(n÷2+number_points_from_mid))
    rho = density(g)
    M = estimated_holder_constant(g, rule, mid_points, rho)
    return select_number_node_per_block(g, OracleM(M))
end

function estimated_holder_constant(g::Observations, ::EstimatedEigenvalues, points, rho)
    eig_res = eigs(g.graph, nev = 1, which = :LM)
    u = eig_res.vectors
    mult = eig_res.values[1]
    return _approx_m_from_delta_f(u, mult, points, rho)
end

function estimated_holder_constant(g::Observations, ::EstimatedDegrees, points, rho)
    d = degree(g.graph)
    mult = (d' * g.graph * d) / (sum(d .^ 2))^2
    return _approx_m_from_delta_f(d, mult, points, rho)
end


function _approx_m_from_delta_f(u, mult, midpoints, ρ, α=1.0)
    sort!(u,dims=1)
    uMid = u[midpoints]
    β₀, β₁  = hcat(ones(length(uMid)), 1:length(uMid)) \ uMid
    h = 2^(α+1) * α * mult^2 * (β₁ * length(uMid)/2 + β₀)^2 * β₁^2 * ρ^(-1/(2*(α+1)))
    return h^(-1/(2*(α+1)))
end

"""
    oracle_bandwidth(A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))

Oracle bandwidth selection for graph histogram, using

```math
\\widehat{h^*}=\\left(2\\left(\\left(d^T d\\right)^{+}\\right)^2 d^T A d \\cdot \\hat{m}
\\hat{b}\\right)^{-\\frac{1}{2}} \\hat{\\rho}_n^{\\frac{1}{4}},
```

where ``d`` is the vector of degree sorted in increasing order,``\\hat{\\rho}_n`` is the
empirical edge density, and  ``m``, ``b`` are the slope and intercept fitted on
``d[n/2-c\\sqrt{n}:n/2+c\\sqrt{n}]`` for some ``c``.
"""
function oracle_bandwidth(
        A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))
    if type ∉ ["eigs", "degs"]
        error("Invalid input type $(type)")
    end

    if alpha != 1
        error("Currently only supports alpha = 1")
    end

    n = size(A, 1)
    midPt = collect(max(
        1, round(Int, (n÷2-c*sqrt(n)))):round(Int,
        (n÷2+c*sqrt(n))))
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
         (lmfit_coef[2] * length(uMid) / 2 + lmfit_coef[1])^2 *
         lmfit_coef[2]^2 *
         rhoHat_inv)^(-1 / (2 * (alpha + 1)))
    return h[1]
end
