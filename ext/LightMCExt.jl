module LightMCExt

using StaticArrays
using Accessors
using NetworkHistogram
import NetworkHistogram: SuffStats, add_sample, remove_sample, make_k_block, score,
                         to_params, AbstractConvertor, to_distribution, get_convertor

using LightMC: DiscreteMarkovChain, SampleChain, transition_matrix, ConvertBinaryMC

# need to define a convertor that only look at the possible transitions and not all of them
struct McConvertor <: AbstractConvertor end

get_convertor(::Val{:mc}; kwargs...) = McConvertor()

function (c::McConvertor)(chain::SampleChain)
    return SVector([SVector(c...) for c in eachcol(chain.transitions)]...)
end

function to_distribution(::McConvertor, transition_matrix; kwargs...)
    return DiscreteMarkovChain(transition_matrix, sum(transition_matrix; dims = 2))
end
struct McSuffStats{M, T} <: SuffStats
    h::SVector{M, T}
end

function McSuffStats(num_states::Int)
    inter = @SVector zeros(SVector{num_states, Int}, num_states)
    return McSuffStats(inter)
end

function add_sample(ss::McSuffStats, sample)
    @inbounds for (i, s) in enumerate(sample)
        ss = @set ss.h[i] = ss.h[i] + s
    end
    return ss
end

function remove_sample(ss::McSuffStats, sample)
    @inbounds for (i, s) in enumerate(sample)
        ss = @set ss.h[i] = ss.h[i] - s
    end
    return ss
end

function _score(counts::SVector)
    n = sum(counts)
    norm_ = max(n, 1)
    return (n - sum(abs2, counts) / norm_) / norm_
end

function score(ss::McSuffStats)
    return sum(_score, ss.h)
end

function to_params(ss::McSuffStats)
    return reduce(hcat, ss.h)
end

end
