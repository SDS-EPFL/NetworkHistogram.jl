include("SymArray.jl")
include("config_rules/include.jl")

function ordered_start_labels(n::Int, k::Int)
    labels = Vector{Int}(undef, n)
    base_size = n ÷ k
    remainder = n % k
    for group in 1:k
        fill!(view(labels, ((group - 1) * base_size + 1):(group * base_size)), group)
    end
    if remainder > 0
        fill!(view(labels, (k * base_size + 1):(k * base_size + remainder)), k)
    end
    return labels
end
