abstract type KSelectionRule end
struct OracleK <: KSelectionRule
    K::Int
end
struct OracleM{F} <: KSelectionRule
    M::F
    α::F
end

function OracleM(M)
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
    if rule.K > number_nodes(g) ÷ 2
        throw(ArgumentError("The number of blocks $(rule.K) is too large for the number \
        of nodes $(number_nodes(g)), it should be at most $(number_nodes(g)÷2)"))
    end
    return rule.K
end

function select_number_node_per_block(g::Observations, rule::OracleM)
    rho = density(g)
    n = number_nodes(g)
    k = max(2, round(Int, (2 * rule.M * rho)^(-1 / 4) * sqrt(n)))
    return select_number_node_per_block(g, OracleK(k))
end

function select_number_node_per_block(g::Observations, rule::EstimatedM)
    n = number_nodes(g)
    c = min(4, sqrt(n) / 8)
    number_points_from_mid = round(Int, c * sqrt(n))
    mid_points = max(1, n ÷ 2 - number_points_from_mid):(n ÷ 2 + number_points_from_mid)
    m = estimated_number_nodes_per_block(g, rule, mid_points, density(g))
    return select_number_node_per_block(g, OracleK(m))
end

function estimated_number_nodes_per_block(
        g::Observations, ::EstimatedEigenvalues, points, rho)
    λ, u = Arpack.eigs(get_adj(g), nev = 1, which = :LM)
    return _approx_k_from_delta_f(u, λ[1], points, rho)
end

function estimated_number_nodes_per_block(
        g::Observations, ::EstimatedDegrees, points, rho)
    d = get_degree(g)
    mult = ((d' * get_adj(g) * d) / (sum(d .^ 2))^2)[1]
    return _approx_k_from_delta_f(d, mult, points, rho)
end

function _approx_k_from_delta_f(u, mult, midpoints, ρ, α = 1.0)
    sort!(u, dims = 1)
    uMid = u[midpoints]
    β₀, β₁ = hcat(ones(length(uMid)), 1:length(uMid)) \ uMid
    # from Olhede and Wolfe (2014), equation (11)
    h = (2^(α + 1) * α * mult^2 * (β₁ * length(uMid) / 2 + β₀)^2 * β₁^2 *
         ρ^(-1))^(-1 / (2 * (α + 1)))
    return round(Int, h)
end
