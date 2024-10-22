import MLJModelInterface
const MMI = MLJModelInterface
using PermutationSymmetricTensors
using Distributions

MMI.@mlj_model mutable struct SBM <: MMI.Probabilistic
    k::Int = 1::(_ > 0)
    D::Val{<:Distribution} = Val{Bernoulli}()
end

function MMI.fit(model::SBM, X, y)
    return model
end
