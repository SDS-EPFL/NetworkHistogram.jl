mutable struct CategoricalSwap{M,F} <: Swap
    index1::Int
    index2::Int
    realized::Matrix{MVector{M, Int}}
    estimated_theta::Matrix{MVector{M, F}}
    log_likelihood::F
end
