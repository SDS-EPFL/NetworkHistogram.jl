# if S is Int, assume the states are ordered and sequential
struct DiscreteMarkovChain{S, T}
    states::Vector{S}
    transitions::Matrix{T}
end

struct SampleChain{S}
    states::Vector{S}
    indices::Vector{Int}
    transitions::Matrix{Int}
end

function state_index(mc::DiscreteMarkovChain{S, T}, state::S) where {S, T}
    findfirst(isequal(state), mc.states)
end

state_space(mc::DiscreteMarkovChain) = mc.states
transition_matrix(mc::DiscreteMarkovChain) = mc.transitions

function stationary_dist(mc::DiscreteMarkovChain)
    T = transition_matrix(mc)
    n = length(state_space(mc))
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

function sample_indices(mc::DiscreteMarkovChain{S, T}, t::Int) where {S, T}
    indices = Vector{Int}(undef, t)
    indices[1] = rand(Categorical(stationary_dist(mc)))
    tr_transposed = transpose(mc.transitions)
    for i in 2:t
        indices[i] = rand(Categorical(tr_transposed[:, indices[i - 1]]))
    end
    return indices
end

function sample(mc::DiscreteMarkovChain{S, T}, t::Int) where {S, T}
    indices = sample_indices(mc, t)
    states = mc.states[indices]
    counts = zeros(Int, length(mc.states), length(mc.states))
    for i in 1:(length(indices) - 1)
        counts[indices[i], indices[i + 1]] += 1
    end
    return SampleChain(states, indices, counts)
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

function loglikelihood(mc::DiscreteMarkovChain{S, T}, chain::Vector{Int}) where {S, T}
    Tr = transition_matrix(mc)
    loglik = log(stationary_dist(mc)[chain[1]])
    for i in 1:(length(chain) - 1)
        loglik += log(Tr[chain[i], chain[i + 1]])
    end
    return loglik
end

function loglikelihood(mc::DiscreteMarkovChain{S, T}, chain::Vector{S}) where {S, T}
    return loglikelihood(mc, state_index.(Ref(mc), chain))
end

function loglikelihood(mc::DiscreteMarkovChain{S, T}, chain::SampleChain{S}) where {S, T}
    return sum(xlogy.(chain.transitions, mc.transitions)) +
           log(stationary_dist(mc)[chain.indices[1]])
end


# user responsability to have the same states...
function fit(mc::DiscreteMarkovChain{S, T}, chain::SampleChain{S}) where {S, T}
    return DiscreteMarkovChain(mc.states, chain.transitions ./ sum(chain.transitions, dims=2))
end
