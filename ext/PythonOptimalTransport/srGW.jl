# Julia implementation of semi-relaxed Gromov-Wasserstein algorithms
# Converted from Python implementation by cvincentcuaz
# This module provides conditional gradient, mirror descent, and MM algorithms
# for semi-relaxed (fused) Gromov-Wasserstein optimal transport

using LinearAlgebra
using Random

# =============================================================================
# Utility Functions
# =============================================================================

# Initialize transport plan for semi-relaxed GW
# Arguments:
#   init_mode: "product", "random", or "random_product"
#   p: source distribution (N1,)
#   N1, N2: dimensions
#   seed: random seed (nothing for no seeding)
# Returns: T - initial transport plan (N1, N2)
function initializer_semirelaxed_GW(
        init_mode::String, p::AbstractVector{T}, N1::Int, N2::Int;
        seed::Union{Int, Nothing} = 0) where {T <: Real}
    if init_mode == "product"
        q = ones(T, N2) / N2
        T_plan = p * q'
    elseif init_mode == "random"
        if !isnothing(seed)
            Random.seed!(seed)
        end
        T_plan = rand(T, N1, N2)
        # Scale to satisfy first marginal constraint
        scale = p ./ sum(T_plan, dims = 2)
        T_plan .*= scale
    elseif init_mode == "random_product"
        if !isnothing(seed)
            Random.seed!(seed)
        end
        q = rand(T, N2)
        q ./= sum(q)
        T_plan = p * q'
    else
        error("Unknown init mode: $init_mode")
    end
    return T_plan
end

# Initialize matrices for symmetric GW computation
function init_matrix_GW2(C1::AbstractMatrix{T}, C2::AbstractMatrix{T},
        p::AbstractVector{T}, q::AbstractVector{T},
        ones_p::AbstractVector{T}, ones_q::AbstractVector{T}) where {T <: Real}
    f1_ = C1 .^ 2
    f2_ = C2 .^ 2
    constC1 = f1_ * (p * ones_q')
    constC2 = (ones_p * q') * f2_
    constC = constC1 + constC2
    hC1 = C1
    hC2 = 2 * C2
    return constC, hC1, hC2
end

# Initialize matrices for asymmetric GW computation
function init_matrix_asymGW2(C1::AbstractMatrix{T}, C2::AbstractMatrix{T},
        p::AbstractVector{T}, q::AbstractVector{T},
        ones_p::AbstractVector{T}, ones_q::AbstractVector{T}) where {T <: Real}
    f1_ = (C1 .^ 2) / 2.0
    f2_ = (C2 .^ 2) / 2.0
    constC1 = f1_ * (p * ones_q')
    constC2 = (ones_p * q') * f2_'
    constC = constC1 + constC2
    hC1 = C1
    hC2 = C2
    return constC, hC1, hC2
end

# Compute tensor product for GW distance
function tensor_product(constC::AbstractMatrix{T}, hC1::AbstractMatrix{T},
        hC2::AbstractMatrix{T}, T_plan::AbstractMatrix{T}) where {T <: Real}
    A = -hC1 * T_plan * hC2'
    return constC + A
end

# =============================================================================
# Conditional Gradient Descent Algorithms
# =============================================================================

# Conditional gradient algorithm for semi-relaxed (fused) Gromov-Wasserstein
# Solves: min_T  α * ⟨L(C₁, C₂) ⊗ T, T⟩ + ⟨M, T⟩
function cg_semirelaxed(C1::AbstractMatrix{T}, p::AbstractVector{T}, C2::AbstractMatrix{T};
        alpha::Real = 1.0, linear_cost::Union{Nothing, AbstractMatrix{T}} = nothing,
        init_mode::String = "product", T_init::Union{Nothing, AbstractMatrix{T}} = nothing,
        symmetry::Bool = true, use_log::Bool = false, eps::Real = 1e-5,
        max_iter::Int = 1000, seed::Int = 0, verbose::Bool = false) where {T <: Real}
    N1, N2 = size(C1, 1), size(C2, 1)

    # Initialize transport plan
    if isnothing(T_init)
        T_plan = initializer_semirelaxed_GW(init_mode, p, N1, N2; seed = seed)
    else
        @assert size(T_init) == (N1, N2)
        T_plan = copy(T_init)
    end

    # Check symmetry
    if isnothing(symmetry)
        symmetry = (C1 == C1') && (C2 == C2')
    end

    # Initialize
    q = vec(sum(T_plan, dims = 1))
    ones_p = ones(T, N1)
    ones_q = ones(T, N2)

    # Compute initial gradient
    if symmetry
        constC, hC1, hC2 = init_matrix_GW2(C1, C2, p, q, ones_p, ones_q)
        G = 2 * tensor_product(constC, hC1, hC2, T_plan)
    else
        constC, hC1, hC2 = init_matrix_asymGW2(C1, C2, p, q, ones_p, ones_q)
        constCt, hC1t, hC2t = init_matrix_asymGW2(C1', C2', p, q, ones_p, ones_q)
        subG = tensor_product(constC, hC1, hC2, T_plan)
        subGt = tensor_product(constCt, hC1t, hC2t, T_plan)
        G = subG + subGt
    end
    G .*= alpha

    srgw_loss = 0.5 * sum(G .* T_plan)
    add_linear_cost = !isnothing(linear_cost)

    if add_linear_cost
        linear_loss = sum(linear_cost .* T_plan)
        current_loss = srgw_loss + linear_loss
        G .+= linear_cost
    else
        current_loss = srgw_loss
    end

    log = use_log ? Dict("loss" => [current_loss]) : nothing
    convergence_criterion = Inf
    outer_count = 0

    while convergence_criterion > eps && outer_count < max_iter
        previous_loss = current_loss

        # Direction finding by solving subproblem on rows
        min_vals = minimum(G, dims = 2)
        X = (G .== min_vals) .* T(1.0)
        row_sums = vec(sum(X, dims = 2))
        scale = p ./ row_sums
        X .*= scale

        # Exact line search
        qX = vec(sum(X, dims = 1))

        if symmetry
            constCX, hC1X, hC2X = init_matrix_GW2(C1, C2, p, qX, ones_p, ones_q)
            GX = 2 * alpha * tensor_product(constCX, hC1X, hC2X, X)
            GXX = 0.5 * sum(GX .* X)
            GXT = 0.5 * sum(GX .* T_plan)

            a = srgw_loss + GXX - 2 * GXT
            b = 2 * (GXT - srgw_loss)
        else
            constCX, hC1X, hC2X = init_matrix_asymGW2(C1, C2, p, qX, ones_p, ones_q)
            constCXt, hC1Xt, hC2Xt = init_matrix_asymGW2(C1', C2', p, qX, ones_p, ones_q)
            subGX = tensor_product(constCX, hC1X, hC2X, X)
            subGXt = tensor_product(constCXt, hC1Xt, hC2Xt, X)
            GX = alpha * (subGX + subGXt)
            GXX = 0.5 * sum(GX .* X)
            subGXt_dotT = sum(subGXt .* T_plan)
            subGTt_dotX = sum(subGt .* X)

            a = srgw_loss + GXX - subGXt_dotT - subGTt_dotX
            b = -2 * srgw_loss + subGXt_dotT + subGTt_dotX
        end

        if add_linear_cost
            linear_loss_X = sum(linear_cost .* X)
            b += linear_loss_X - linear_loss
        end

        # Compute step size
        if a > 0
            gamma = min(1, max(0, -b / (2 * a)))
        elseif a + b < 0
            gamma = 1
        else
            gamma = 0
        end

        # Update
        T_plan .= (1 - gamma) * T_plan + gamma * X
        current_loss += a * gamma^2 + b * gamma

        if add_linear_cost
            linear_loss = (1 - gamma) * linear_loss + gamma * linear_loss_X
            srgw_loss = current_loss - linear_loss
            G .= (1 - gamma) * G + gamma * (GX + linear_cost)
        else
            srgw_loss = current_loss
            G .= (1 - gamma) * G + gamma * GX
        end

        outer_count += 1
        use_log && push!(log["loss"], current_loss)

        convergence_criterion = abs(previous_loss - current_loss) /
                                (abs(previous_loss) + 1e-15)
    end

    return use_log ? (T_plan, current_loss, log) : (T_plan, current_loss)
end

# Conditional gradient for semi-relaxed Gromov-Wasserstein
# Wrapper for cg_semirelaxed with α=1 and no linear cost
function cg_semirelaxed_gromov_wasserstein(C1::AbstractMatrix{T}, p::AbstractVector{T},
        C2::AbstractMatrix{T}; kwargs...) where {T <: Real}
    return cg_semirelaxed(C1, p, C2; alpha = 1.0, linear_cost = nothing, kwargs...)
end

# Conditional gradient for semi-relaxed fused Gromov-Wasserstein
# A1, A2: Feature matrices (N1×d), (N2×d)
# alpha: Trade-off parameter (0 for pure OT, 1 for pure GW)
function cg_semirelaxed_fused_gromov_wasserstein(
        C1::AbstractMatrix{T}, A1::AbstractMatrix{T},
        p::AbstractVector{T}, C2::AbstractMatrix{T},
        A2::AbstractMatrix{T}, alpha::Real;
        kwargs...) where {T <: Real}
    N1, N2 = size(A1, 1), size(A2, 1)
    d = size(A1, 2)

    # Compute Euclidean distance matrix between features
    A1_sq = sum(A1 .^ 2, dims = 2) * ones(T, 1, N2)
    A2_sq = ones(T, N1, 1) * sum(A2 .^ 2, dims = 2)'
    D = A1_sq + A2_sq - 2 * A1 * A2'

    return cg_semirelaxed(
        C1, p, C2; alpha = alpha, linear_cost = (1 - alpha) * D, kwargs...)
end

# =============================================================================
# Mirror Descent Algorithms (Entropic Regularization)
# =============================================================================

# Mirror descent algorithm using KL geometry for semi-relaxed (fused) GW
# gamma_entropy: Entropy regularization parameter (must be > 0)
function md_semirelaxed(C1::AbstractMatrix{T}, p::AbstractVector{T}, C2::AbstractMatrix{T},
        gamma_entropy::Real; alpha::Real = 1.0,
        linear_cost::Union{Nothing, AbstractMatrix{T}} = nothing,
        init_mode::String = "product", T_init::Union{Nothing, AbstractMatrix{T}} = nothing,
        symmetry::Bool = true, use_log::Bool = false, eps::Real = 1e-5,
        max_iter::Int = 1000, seed::Int = 0, verbose::Bool = false) where {T <: Real}
    @assert gamma_entropy>0 "gamma_entropy must be positive"

    N1, N2 = size(C1, 1), size(C2, 1)

    # Initialize transport plan
    if isnothing(T_init)
        T_plan = initializer_semirelaxed_GW(init_mode, p, N1, N2; seed = seed)
    else
        @assert size(T_init) == (N1, N2)
        T_plan = copy(T_init)
    end

    # Check symmetry
    if isnothing(symmetry)
        symmetry = (C1 == C1') && (C2 == C2')
    end

    # Initialize
    q = vec(sum(T_plan, dims = 1))
    ones_p = ones(T, N1)
    ones_q = ones(T, N2)

    # Compute initial gradient
    if symmetry
        constC, hC1, hC2 = init_matrix_GW2(C1, C2, p, q, ones_p, ones_q)
        G = 2 * alpha * tensor_product(constC, hC1, hC2, T_plan)
    else
        constC, hC1, hC2 = init_matrix_asymGW2(C1, C2, p, q, ones_p, ones_q)
        constCt, hC1t, hC2t = init_matrix_asymGW2(C1', C2', p, q, ones_p, ones_q)
        subG = tensor_product(constC, hC1, hC2, T_plan)
        subGt = tensor_product(constCt, hC1t, hC2t, T_plan)
        G = alpha * (subG + subGt)
    end

    current_loss = 0.5 * sum(G .* T_plan)
    add_linear_cost = !isnothing(linear_cost)

    if add_linear_cost
        linear_loss = sum(linear_cost .* T_plan)
        current_loss += linear_loss
        G .+= linear_cost
    end

    log = use_log ? Dict("loss" => [current_loss]) : nothing
    convergence_criterion = Inf
    outer_count = 0

    while convergence_criterion > eps && outer_count < max_iter
        previous_loss = current_loss

        # Compute Bregman projection
        M = G - gamma_entropy * Base.log.(T_plan)
        K = Base.exp.(-M / gamma_entropy)
        scaling = p ./ vec(sum(K, dims = 2))
        T_plan .= (scaling .* ones(T, 1, N2)) .* K

        q = vec(sum(T_plan, dims = 1))

        # Update gradient
        if symmetry
            constC, hC1, hC2 = init_matrix_GW2(C1, C2, p, q, ones_p, ones_q)
            G = 2 * alpha * tensor_product(constC, hC1, hC2, T_plan)
        else
            constC, hC1, hC2 = init_matrix_asymGW2(C1, C2, p, q, ones_p, ones_q)
            constCt, hC1t, hC2t = init_matrix_asymGW2(C1', C2', p, q, ones_p, ones_q)
            subG = tensor_product(constC, hC1, hC2, T_plan)
            subGt = tensor_product(constCt, hC1t, hC2t, T_plan)
            G = alpha * (subG + subGt)
        end

        current_loss = 0.5 * sum(G .* T_plan)

        if add_linear_cost
            linear_loss = sum(linear_cost .* T_plan)
            current_loss += linear_loss
            G .+= linear_cost
        end

        outer_count += 1
        use_log && push!(log["loss"], current_loss)

        convergence_criterion = abs(previous_loss - current_loss) /
                                (abs(previous_loss) + 1e-15)
    end

    return use_log ? (T_plan, current_loss, log) : (T_plan, current_loss)
end

# Mirror descent for semi-relaxed Gromov-Wasserstein with entropic regularization
function md_semirelaxed_gromov_wasserstein(C1::AbstractMatrix{T}, p::AbstractVector{T},
        C2::AbstractMatrix{T}, gamma_entropy::Real;
        kwargs...) where {T <: Real}
    return md_semirelaxed(
        C1, p, C2, gamma_entropy; alpha = 1.0, linear_cost = nothing, kwargs...)
end

# Mirror descent for semi-relaxed fused Gromov-Wasserstein with entropic regularization
function md_semirelaxed_fused_gromov_wasserstein(
        C1::AbstractMatrix{T}, A1::AbstractMatrix{T},
        p::AbstractVector{T}, C2::AbstractMatrix{T},
        A2::AbstractMatrix{T}, gamma_entropy::Real,
        alpha::Real; kwargs...) where {T <: Real}
    N1, N2 = size(A1, 1), size(A2, 1)
    d = size(A1, 2)

    # Compute Euclidean distance matrix
    A1_sq = sum(A1 .^ 2, dims = 2) * ones(T, 1, N2)
    A2_sq = ones(T, N1, 1) * sum(A2 .^ 2, dims = 2)'
    D = A1_sq + A2_sq - 2 * A1 * A2'

    return md_semirelaxed(C1, p, C2, gamma_entropy; alpha = alpha,
        linear_cost = (1 - alpha) * D, kwargs...)
end

# =============================================================================
# Majorization-Minimization Algorithms with Sparsity Regularization
# =============================================================================

# MM algorithm with ℓₚ-ℓ₁ sparsity regularization for semi-relaxed (fused) GW
# Solves: min_T  α⟨L(C₁,C₂)⊗T,T⟩ + ⟨M,T⟩ + λ∑ⱼ(∑ᵢTᵢⱼ)^p
function mm_lpl1_semirelaxed(
        C1::AbstractMatrix{T}, p::AbstractVector{T}, C2::AbstractMatrix{T},
        gamma_entropy::Real; alpha::Real = 1.0,
        linear_cost::Union{Nothing, AbstractMatrix{T}} = nothing,
        T_init::Union{Nothing, AbstractMatrix{T}} = nothing,
        init_mode::String = "product", symmetry::Bool = true,
        p_reg::Real = 0.5, lambda_reg::Real = 0.001,
        use_log::Bool = false, use_warmstart::Bool = false,
        eps_inner::Real = 1e-6, eps_outer::Real = 1e-6,
        max_iter_inner::Int = 1000, max_iter_outer::Int = 50,
        seed::Int = 0, verbose::Bool = false,
        inner_log::Bool = false) where {T <: Real}
    @assert 0<p_reg<1 "p_reg must be in (0, 1)"
    @assert gamma_entropy>=0 "gamma_entropy must be non-negative"

    N1, N2 = size(C1, 1), size(C2, 1)

    # Initialize
    if isnothing(T_init)
        T_plan = initializer_semirelaxed_GW(init_mode, p, N1, N2; seed = seed)
        T_init_warm = use_warmstart ? copy(T_plan) : nothing
    else
        @assert size(T_init) == (N1, N2)
        T_plan = copy(T_init)
        T_init_warm = nothing
    end

    # Inner solver selection
    if gamma_entropy == 0
        inner_solver = (total_linear_cost,
            T_init_local) -> cg_semirelaxed(
            C1, p, C2; alpha = alpha, linear_cost = total_linear_cost,
            init_mode = init_mode, T_init = T_init_local, symmetry = symmetry,
            use_log = inner_log, eps = eps_inner, max_iter = max_iter_inner,
            seed = seed, verbose = verbose
        )
    else
        inner_solver = (total_linear_cost,
            T_init_local) -> md_semirelaxed(
            C1, p, C2, gamma_entropy; alpha = alpha, linear_cost = total_linear_cost,
            init_mode = init_mode, T_init = T_init_local, symmetry = symmetry,
            use_log = inner_log, eps = eps_inner, max_iter = max_iter_inner,
            seed = seed, verbose = verbose
        )
    end

    reg_linear_cost = zeros(T, N1, N2)
    total_linear_cost = isnothing(linear_cost) ? nothing : copy(linear_cost)

    best_T = copy(T_plan)
    ones_p = ones(T, N1, 1)

    log = use_log ? Dict("loss" => T[], "inner_loss" => []) : nothing
    best_loss = T(Inf)
    current_loss = T(1e15)
    convergence_criterion = Inf
    outer_count = 0

    while convergence_criterion > eps_outer && outer_count < max_iter_outer
        previous_loss = current_loss

        # Solve generalized problem
        result = inner_solver(total_linear_cost, use_warmstart ? T_init_warm : nothing)

        if inner_log
            T_plan, majorization_loss, inner_log_data = result
        else
            T_plan, majorization_loss = result
        end

        # Compute linearized reg loss
        linearized_reg_loss = sum(reg_linear_cost .* T_plan)

        if use_warmstart
            T_init_warm = copy(T_plan)
        end

        # Update regularization
        q = vec(sum(T_plan, dims = 1))
        reg_loss = lambda_reg * sum((q .+ 1e-15) .^ p_reg)
        current_loss = majorization_loss - linearized_reg_loss + reg_loss

        reg_linear_cost .= lambda_reg * p_reg * ((ones_p * q') .+ 1e-15) .^ (p_reg - 1.0)

        if isnothing(linear_cost)
            total_linear_cost = reg_linear_cost
        else
            total_linear_cost = reg_linear_cost + linear_cost
        end

        if verbose
            println("Outer iter $outer_count: loss = $current_loss, q = $q")
        end

        outer_count += 1

        if use_log
            push!(log["loss"], current_loss)
            inner_log && push!(log["inner_loss"], inner_log_data)
        end

        convergence_criterion = abs(previous_loss - current_loss) /
                                (abs(previous_loss) + 1e-15)

        if current_loss < best_loss
            best_loss = current_loss
            best_T = copy(T_plan)
        end
    end

    return use_log ? (best_T, best_loss, log) : (best_T, best_loss)
end

# MM algorithm with sparsity for semi-relaxed Gromov-Wasserstein
function mm_lpl1_semirelaxed_gromov_wasserstein(
        C1::AbstractMatrix{T}, p::AbstractVector{T},
        C2::AbstractMatrix{T}, gamma_entropy::Real;
        kwargs...) where {T <: Real}
    return mm_lpl1_semirelaxed(C1, p, C2, gamma_entropy; alpha = 1.0,
        linear_cost = nothing, kwargs...)
end

# MM algorithm with sparsity for semi-relaxed fused Gromov-Wasserstein
function mm_lpl1_semirelaxed_fused_gromov_wasserstein(
        C1::AbstractMatrix{T}, A1::AbstractMatrix{T},
        p::AbstractVector{T}, C2::AbstractMatrix{T},
        A2::AbstractMatrix{T}, alpha::Real,
        gamma_entropy::Real; kwargs...) where {T <: Real}
    N1, N2 = size(A1, 1), size(A2, 1)
    d = size(A1, 2)

    # Compute Euclidean distance matrix
    A1_sq = sum(A1 .^ 2, dims = 2) * ones(T, 1, N2)
    A2_sq = ones(T, N1, 1) * sum(A2 .^ 2, dims = 2)'
    D = A1_sq + A2_sq - 2 * A1 * A2'

    return mm_lpl1_semirelaxed(C1, p, C2, gamma_entropy; alpha = alpha,
        linear_cost = (1 - alpha) * D, kwargs...)
end
