module LightMCExt

using NetworkHistogram
using LightMC

import NetworkHistogram: agg_params, logpdf, sample, params, distance, _fast_compressed_obs,
                         from_adjs_to_decorated
using LightMC: DiscreteMarkovChain, SampleChain, transition_matrix, ConvertBinaryMC

logpdf(d::DiscreteMarkovChain, x) = LightMC.logpdf(d, x)
sample(x::DiscreteMarkovChain, args...) = LightMC.sample(x, args...)
params(d) = LightMC.params(d)

function agg_params(d1::DiscreteMarkovChain, d2::DiscreteMarkovChain, w1, w2)
    s1 = Int(sign(w1))
    s2 = Int(sign(w2))
    return DiscreteMarkovChain(s1 .* d1.transitions .+ s2 .* d2.transitions,
        s1 .* d1.normalization .+ s2 .* d2.normalization)
end

function distance(d1::DiscreteMarkovChain, d2::DiscreteMarkovChain)
    mean(x -> x^2, transition_matrix(d1) - transition_matrix(d2))
end
function distance(d1::SampleChain, d2::SampleChain)
    mean(x -> x^2, transition_matrix(d1) - transition_matrix(d2))
end
params(d::DiscreteMarkovChain) = (d.transitions, d.normalization)

function _fast_compressed_obs(d::DiscreteMarkovChain, x::SampleChain, zeroinflated)
    return x
end

function from_adjs_to_decorated(adjs::AbstractArray{T, 3}, converter::ConvertBinaryMC,
        threshold = 0.0) where {T <: Union{Missing, Real}}
    sample_chain = MC.periodic_chain(adjs[1, 4, :], converter)
    graph = Matrix{Union{typeof(sample_chain), Missing}}(
        undef, size(adjs, 1), size(adjs, 2))
    counts_t = sum(adjs, dims = 3)
    for j in axes(adjs, 2)
        for i in axes(adjs, 1)
            if i == j || counts_t[i, j] <= threshold * size(adjs, 3)
                graph[i, j] = missing
            else
                graph[i, j] = LightMC.periodic_chain(adjs[i, j, :], converter)
            end
        end
    end
    return graph
end

function from_adjs_to_decorated(adjs::AbstractArray{T, 2}, converter::ConvertBinaryMC,
        threshold = 0.0) where {T <: Union{Missing, AbstractArray}}
    sample_chain = LightMC.periodic_chain(adjs[1, 4], converter)
    graph = Matrix{Union{typeof(sample_chain), Missing}}(
        undef, size(adjs, 1), size(adjs, 2))
    for j in axes(adjs, 2)
        for i in axes(adjs, 1)
            if i == j || sum(adjs[i, j]) <= threshold * length(adjs[i, j])
                graph[i, j] = missing
            else
                graph[i, j] = LightMC.periodic_chain(adjs[i, j], converter)
            end
        end
    end
    return graph
end

end
