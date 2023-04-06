function checkadjacency(A)
    @assert eltype(A) <: Real
    if !(eltype(A) === Bool)
        @assert all(a ∈ [zero(eltype(A)), one(eltype(A))] for a in A) "All elements of the ajacency matrix should be zero or one."
    end
    @assert issymmetric(A)
    @assert all(A[i, i] == zero(eltype(A)) for i in 1:size(A, 1)) "The diagonal of the adjacency matrix should all be zeros."
    return nothing
end

"""
    graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule = RandomNodeSwap(),
    starting_assignment_rule = RandomStart(), accept_rule = Strict(),
    stop_rule = PreviousBestValue(3), record_trace=true)

Compute the graph histogram.

# Arguments
TBW
"""
function graphhist(A; h = select_bandwidth(A), maxitr = 1000,
                   swap_rule::NodeSwapRule = RandomNodeSwap(),
                   starting_assignment_rule::StartingAssignment = RandomStart(),
                   accept_rule::AcceptRule = Strict(),
                   stop_rule::StopRule = PreviousBestValue(3), record_trace = true)
    checkadjacency(A)
    @assert maxitr > 0

    return _graphhist(A, Val{record_trace}(), h = h, maxitr = maxitr, swap_rule = swap_rule,
                      starting_assignment_rule = starting_assignment_rule,
                      accept_rule = accept_rule,
                      stop_rule = stop_rule)
end

"""
    _graphhist(A, record_trace=Val{true}(); h, maxitr, swap_rule, starting_assignment_rule, accept_rule, stop_rule)

Internal version of `graphhist` which is type stable.
"""
function _graphhist(A, record_trace = Val{true}(); h, maxitr, swap_rule,
                    starting_assignment_rule, accept_rule, stop_rule)
    best, current, proposal, history = initialize(A, h, starting_assignment_rule,
                                                  record_trace)

    for i in 1:maxitr
        proposal = create_proposal!(history, i, proposal, current, A, swap_rule)
        current = accept_reject_update!(history, i, proposal, current, accept_rule)
        best = update_best!(history, i, current, best)
        if stopping_rule(history, stop_rule)
            break
        end
    end

    return graphhist_format_output(best, history)
end

"""
    graphhist_format_output(best, history)

Formates the `graphhist` output depending on the type of `history` requested by the user.
"""
function graphhist_format_output(best, history::TraceHistory)
    return (graphhist = GraphHist(best), trace = history, likelihood = best.likelihood)
end
function graphhist_format_output(best, history::NoTraceHistory)
    return (graphhist = GraphHist(best), likelihood = history.best_likelihood)
end

function update_best!(history::GraphOptimizationHistory, iteration::Int,
                      current::Assignment,
                      best::Assignment)
    if current.likelihood > best.likelihood
        update_best!(history, iteration, current.likelihood)
        deepcopy!(best, current)
    end
    return best
end

"""
    initialize(A, h, starting_assignment_rule, record_trace)

Initialize the memory required for finding optimal graph histogram.
"""
function initialize(A, h, starting_assignment_rule, record_trace)
    node_labels, group_size = initialise_node_labels(A, h, starting_assignment_rule)
    proposal = Assignment(A, node_labels, group_size)
    current = deepcopy(proposal)
    best = deepcopy(proposal)
    history = initialize_history(best, current, proposal, record_trace)
    return best, current, proposal, history
end

function select_bandwith(A, type = "degs", alpha = 1, c = 1)
    h = oracle_bandwidth(A, type, alpha, c)
    return sanitize_bandwidth(h, size(A, 1))
end

function oracle_bandwidth(A, type = "degs", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))
    if type ∉ ["eigs", "degs"]
        error("Invalid input type $(type)")
    end

    if alpha != 1
        error("Currently only supports alpha = 1")
    end

    n = size(A, 1)
    midPt = collect(round(Int, (n ÷ 2 - c * sqrt(n))):round(Int, (n ÷ 2 + c * sqrt(n))))
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

function sanitize_bandwidth(h::Real, n::Int)::Int
    h = max(2, min(n, round(Int, h)))
    lastGroupSize = n % h
    # step down h, to avoid singleton final group
    while lastGroupSize == 1 && h > 2
        h -= 1
        lastGroupSize = n % h
    end
    return h
end
