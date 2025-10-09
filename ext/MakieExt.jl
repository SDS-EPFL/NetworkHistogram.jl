module MakieExt

using NetworkHistogram
using Makie

import NetworkHistogram: get_probability_matrix, Assignment, heatmap_params, params,
                         number_nodes
import Distributions: params

_splatter_args(ps) = vcat(vec.(ps)...)
_extract_params(d) = _splatter_args(params(d))

function Makie.convert_arguments(::Type{<:AbstractPlot}, a::Assignment)
    params_matrix = map(_extract_params, get_probability_matrix(a))
    ps = (getindex.(params_matrix, i) for i in 1:length(params_matrix[1, 2]))
    return ps
end

function heatmap_params(a; colormap = :balance, ordering = false)
    params_matrix = map(_extract_params, get_probability_matrix(a))
    if ordering
        perm = sortperm(a.node_labels)
    else
        perm = 1:number_nodes(a)
    end
    params_matrix = params_matrix[perm, perm]
    num_params = length(params_matrix[1, 2])
    # Compute rows and columns such that rows * columns >= num_params and as square as possible
    rows = floor(Int, sqrt(num_params))
    cols = ceil(Int, num_params / rows)
    if rows * cols < num_params
        rows += 1
    end
    fig = Figure(size = (800, 800))
    # create a grid of subplots with rows x cols cells
    axes = [Axis(fig[i, j]) for i in 1:rows, j in 1:cols]

    for i in 1:num_params
        heatmap!(axes[i], getindex.(params_matrix, i)[perm, perm], colormap = colormap)
        axes[i].title = "Parameter $i"
    end
    return fig
end

export heatmap_params
end
