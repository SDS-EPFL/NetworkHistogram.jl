
"""
    align_partitions(x::AbstractVector{<:Integer}, y::AbstractVector{<:Integer})

Align labels of partition `y` to match partition `x` using optimal matching.

Returns `(y_aligned, mapping)` where:
- `y_aligned`: Vector with same length as `y`, with labels relabeled to match `x`
- `mapping`: Dictionary mapping original `y` labels to aligned `x` labels

The alignment maximizes the overlap between partitions using the Hungarian algorithm.
Unmatched labels from `y` are assigned to unused labels from `x`.

# Arguments
- `x::AbstractVector{<:Integer}`: Reference partition labels
- `y::AbstractVector{<:Integer}`: Partition labels to align

# Examples
```julia
x = [1, 1, 2, 2, 3]
y = [2, 2, 1, 1, 3]
y_aligned, mapping = align_partitions(x, y)
# y_aligned == [1, 1, 2, 2, 3]
# mapping == Dict(2 => 1, 1 => 2, 3 => 3)
```
"""
function align_partitions(x::AbstractVector{<:Integer}, y::AbstractVector{<:Integer})
    @argcheck length(x)==length(y) "Partitions must have same length"

    # Get unique labels and create index mappings
    xlabs = sort!(unique(x))
    ylabs = sort!(unique(y))
    Kx, Ky = length(xlabs), length(ylabs)

    # Build contingency matrix: C[i,j] = count where x==xlabs[i] and y==ylabs[j]
    C = zeros(Int, Kx, Ky)
    for k in eachindex(x, y)
        i = searchsortedfirst(xlabs, x[k])
        j = searchsortedfirst(ylabs, y[k])
        @inbounds C[i, j] = 1
    end

    # Solve maximum weight assignment on padded square matrix
    K = max(Kx, Ky)
    W = zeros(Int, K, K)
    @views W[1:Kx, 1:Ky] .= -C
    assignment, cost = hungarian(W)

    # Build mapping from y to x labels
    mapping = Dict{eltype(ylabs), eltype(xlabs)}()
    used_x = Set{eltype(xlabs)}()

    # Map matched pairs within actual partition sizes
    for i in 1:Kx
        j = assignment[i]
        if 1 ≤ j ≤ Ky
            mapping[ylabs[j]] = xlabs[i]
            push!(used_x, xlabs[i])
        end
    end

    # Assign unmatched y labels to unused x labels
    unused_x = filter(∉(used_x), xlabs)
    for j in 1:Ky
        if !haskey(mapping, ylabs[j])
            @argcheck !isempty(unused_x) "Insufficient x labels for alignment"
            mapping[ylabs[j]] = popfirst!(unused_x)
        end
    end

    # Relabel y using the mapping
    y_aligned = map(ylab -> mapping[ylab], y)

    return y_aligned, mapping
end
