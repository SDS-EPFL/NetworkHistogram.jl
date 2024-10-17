# mutable struct CategoricalData{F}
#     counts::Matrix{Int}
#     realized::Matrix{Int}
#     estimated_theta::Matrix{F}
#     A::Matrix{Int}  # possible improvement by using an adjacency list  Graphs.SimpleGraphs.adj(G)
#     log_likelihood::F
# end

# const CategoricalAssignment{T, F} = Assignment{T, CategoricalData{F}}
# const CategoricalInitRule{S, F} = InitRule{S, Val{CategoricalData}}
