"""
MakieExt - Visualization extension for NetworkHistogram

Provides plotting capabilities for Assignment and BlockModel objects using Makie.jl.
"""
module MakieExt

using NetworkHistogram
using Makie
using StatsBase: countmap

import NetworkHistogram: get_probability_matrix, Assignment, heatmap_params,
                         number_nodes, number_groups, Dist, BlockModel
import Distributions
import StatsAPI

# Helper functions to extract distribution parameters
vec_mine(x) = vec(x)
vec_mine(x::Real) = x

_splatter_args(ps) = vcat(vec_mine.(ps)...)
_extract_params(d) = _splatter_args(StatsAPI.params(d))

"""
    Makie.convert_arguments(::Type{<:AbstractPlot}, a::Assignment)

Convert an Assignment to plottable data by extracting distribution parameters.
"""
function Makie.convert_arguments(::Type{<:AbstractPlot}, a::Assignment)
    params_matrix = map(_extract_params, get_probability_matrix(a))
    ps = (getindex.(params_matrix, i) for i in 1:length(params_matrix[1, 2]))
    return ps
end

"""
    Makie.convert_arguments(::Type{<:Heatmap}, sbm::BlockModel)

Convert a BlockModel to heatmap arguments for Bernoulli distributions.
"""
function Makie.convert_arguments(::Type{<:Heatmap},
        sbm::BlockModel{D}) where {D <: Union{
        Dist{T}, T} where {T <: NetworkHistogram.Bernoulli}}
    return (0:0.01:1, 0:0.01:1, (x, y) -> first(StatsAPI.params(sbm[x, y])))
end

"""
    heatmap_params(a; colormap=:binary, ordering=false, colorrange=nothing, group_match=1:number_groups(a))

Create a heatmap visualization of distribution parameters in an Assignment.

# Arguments
- `a::Assignment`: The assignment to visualize
- `colormap`: Color scheme for the heatmap (default: :binary)
- `ordering::Bool`: Whether to sort nodes by their group labels (default: false)
- `colorrange`: Range for color mapping (default: auto-computed from data)
- `group_match`: Mapping for group indices (default: identity)

# Returns
- `Figure`: Makie figure with heatmap(s) showing distribution parameters

# Note
For multi-parameter distributions, creates a grid of heatmaps, one per parameter.
"""
function heatmap_params(a; colormap = :binary, ordering = false,
        colorrange = nothing, group_match = 1:number_groups(a))
    node_labels_new = map(x -> group_match[x], a.node_labels)

    params_matrix = map(
        _extract_params, get_probability_matrix(a, nothing, node_labels_new))

    if ordering
        perm = sortperm(a.node_labels)
    else
        perm = 1:number_nodes(a)
    end
    params_matrix = params_matrix[perm, perm]

    if isnothing(colorrange)
        colorrange = extrema(_splatter_args(params_matrix))
    end

    num_params = length(params_matrix[1, 2])
    # Compute grid dimensions: make as square as possible
    rows = floor(Int, sqrt(num_params))
    cols = ceil(Int, num_params / rows)
    if rows * cols < num_params
        rows += 1
    end

    default_size = 300
    fig = Figure()

    # Create grid of subplots
    axes = [Axis(fig[i, j], width = default_size, height = default_size)
            for i in 1:rows, j in 1:cols]

    for i in 1:num_params
        heatmap!(axes[i], getindex.(params_matrix, i)[perm, perm],
            colormap = colormap, colorrange = colorrange)
        axes[i].title = "Parameter $i"
    end

    Colorbar(fig[1:rows, cols + 1], limits = colorrange, colormap = colormap,
        label = "Parameter value", width = ceil(Int, sqrt(default_size)))
    resize_to_layout!(fig)
    return fig
end

"""
    order_groups(a::Assignment, latents::AbstractVector)

Order groups based on true latent variables (heuristic alignment).

This is a heuristic approach to match estimated groups to ground truth orderings
by analyzing the overlap between sorted latents and group assignments.

# Arguments
- `a::Assignment`: The assignment with estimated groups
- `latents::AbstractVector`: True latent variables (e.g., block memberships)

# Returns
- Permutation vector for reordering groups
"""
function order_groups(a::Assignment, latents::AbstractVector)
    n = number_nodes(a)
    k = number_groups(a)
    sort_perm = sortperm(latents)
    sorted_group_labels = a.node_labels[sort_perm]
    dummy_group_labels = repeat(1:k, inner = n ÷ k + 1)[1:n]
    counts = Dict(group => countmap(dummy_group_labels[sorted_group_labels .== group])
    for group in 1:k)
    return sort(1:k, by = x -> Tuple(get(counts[x], g, 0) for g in 1:k), rev = true)
end

"""
    align_sbm_true_latents!(sbm::NetworkHistogram.BlockModel, a::Assignment, latents)

Align a BlockModel's groups to match true latent variables.

# Arguments
- `sbm::BlockModel`: The block model to align (modified in-place)
- `a::Assignment`: The assignment
- `latents`: True latent variables

# Note
This modifies `sbm` in-place to reorder its blocks.
"""
function align_sbm_true_latents!(sbm::NetworkHistogram.BlockModel, a::Assignment, latents)
    NetworkHistogram.align_sbm!(sbm, order_groups(a, latents))
end

export heatmap_params
end
