# if S is Int, assume the states are ordered and sequential
# should store everything in transpose, will be faster but way more
# complicated to read
struct DiscreteMarkovChain{S, M <: AbstractMatrix}
    states::Vector{S}
    transitions::M
end

struct SampleChain{S, M <: AbstractMatrix}
    states::Vector{S}
    indices::Vector{Int}
    transitions::M
end

Base.zero(::DiscreteMarkovChain) = DiscreteMarkovChain(Int[], zeros(Int, 0, 0))
Base.zero(::SampleChain{S}) where {S} = SampleChain(S[], Int[], zeros(Int, 1, 1))

function state_index(mc::DiscreteMarkovChain{S}, state::S) where {S}
    findfirst(isequal(state), mc.states)
end

state_space(mc::DiscreteMarkovChain) = mc.states
transition_matrix(mc::DiscreteMarkovChain) = mc.transitions

function stationary_dist(mc::DiscreteMarkovChain)
    T = transition_matrix(mc)
    F = eigen(T')
    tol = 1e-8
    idx = findfirst(abs.(F.values .- 1) .< tol)
    if idx === nothing
        error("No eigenvalue equal (within tolerance) to 1 found. The chain may not be ergodic.")
    end
    # Extract the corresponding eigenvector and normalize it to sum to 1.
    pi = real(F.vectors[:, idx])
    return pi ./ sum(pi)
end

function stationary_dist(mc::DiscreteMarkovChain{S, <:SparseMatrixCSC}) where {S}
    T = transition_matrix(mc)
    vals, vecs, _ = eigsolve(T')
    tol = 1e-8
    idx = findfirst(abs.(vals .- 1) .< tol)
    if idx === nothing
        error("No eigenvalue equal (within tolerance) to 1 found. The chain may not be ergodic.")
    end
    # Extract the corresponding eigenvector and normalize it to sum to 1.
    pi = Real.(vecs[idx])
    result = pi ./ sum(pi)
    return result
end

function sample_indices(mc::DiscreteMarkovChain, t::Int)
    indices = Vector{Int}(undef, t)
    indices[1] = rand(Categorical(stationary_dist(mc)))
    tr_transposed = transpose(mc.transitions)
    for i in 2:t
        indices[i] = rand(Categorical(tr_transposed[:, indices[i - 1]]))
    end
    return indices
end

function sample(mc::DiscreteMarkovChain, t::Int)
    indices = sample_indices(mc, t)
    states = mc.states[indices]
    counts = zeros(Int, length(mc.states), length(mc.states))
    for i in 1:(length(indices) - 1)
        counts[indices[i], indices[i + 1]] += 1
    end
    return SampleChain(states, indices, counts)
end

function sample(mc::DiscreteMarkovChain{S, <:SparseMatrixCSC}, t::Int) where {S}
    indices = sample_indices(mc, t)
    states = mc.states[indices]
    counts = zeros(Int, length(mc.states), length(mc.states))
    for i in 1:(length(indices) - 1)
        counts[indices[i], indices[i + 1]] += 1
    end
    return SampleChain(states, indices, sparse(counts))
end

## yes I know this is awful and does not return a proper chain, but...
function Base.:+(a::DiscreteMarkovChain, b::DiscreteMarkovChain)
    return DiscreteMarkovChain(
        a.states,
        a.transitions .+ b.transitions)
end

function Base.:-(a::DiscreteMarkovChain, b::DiscreteMarkovChain)
    return DiscreteMarkovChain(
        a.states,
        a.transitions .- b.transitions)
end

function Base.:*(a::DiscreteMarkovChain, c::Real)
    return DiscreteMarkovChain(
        a.states,
        a.transitions .* c)
end

Base.:*(c::Real, a::DiscreteMarkovChain) = a * c

function Base.:/(a::DiscreteMarkovChain, c::Real)
    return DiscreteMarkovChain(
        a.states,
        a.transitions ./ c)
end

function loglikelihood(mc::DiscreteMarkovChain{S, M}, chain::Vector{Int}) where {S, M}
    Tr = transition_matrix(mc)
    probas = Vector{Float64}(undef, length(chain))
    probas[1] = stationary_dist(mc)[chain[1]]
    for i in 1:(length(chain) - 1)
        probas[i + 1] = Tr[chain[i], chain[i + 1]]
    end
    return sum(log, probas)
end

function loglikelihood(
        mc::DiscreteMarkovChain{S, M1}, chain::Vector{S}) where {S, M1}
    return loglikelihood(mc, state_index.(Ref(mc), chain))
end

#without the first state, huge computational speedup
function loglikelihood(
        mc::DiscreteMarkovChain{S, M1}, chain::SampleChain{S, M2}) where {S, M1, M2}
    return sum(map(xlogy, chain.transitions, mc.transitions)) #+log(stationary_dist(mc)[chain.indices[1]])
end



# user responsability to have the same states...
function fit(
        mc::DiscreteMarkovChain{S, M1}, chain::SampleChain{S, M2}) where {S, M1, M2}
    return DiscreteMarkovChain(
        mc.states, make_row_stochastic(chain.transitions))
end



function make_row_stochastic(A::M) where {M <: AbstractMatrix}
    f(row) = sum(row) == 0 ? ones(length(row)) / length(row) : row ./ sum(row)
    return mapslices(f, A, dims = 2)
end
