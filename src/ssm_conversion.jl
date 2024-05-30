struct GraphShapeHist{T,M}
    θ::Array{T, M}
    node_labels::Vector{Int}
    edge_labels::Vector{Int}
    num_layers::Int
    num_shapes::Int
    shapes_labels::Vector{Int}
end

function make_x_matrix(block_approx::GraphHist{T,M}) where {T, M}
    k = get_num_blocks(block_approx)
    n = length(block_approx.node_labels)
    max_num_shapes = k * (k + 1) ÷ 2
    X = Matrix{T}(undef, get_num_params_per_block(block_approx), max_num_shapes)
    indices_shape = CartesianIndex.((i, j) for j in 1:k for i in 1:k if i ≥ j)
    for (index, i_cart) in enumerate(indices_shape)
        X[:, index] .= block_approx.θ[i_cart, :]
    end
    indices_edges = CartesianIndex.((i, j) for j in 1:n for i in 1:n if i > j)
    return X, k, n, max_num_shapes, indices_shape, indices_edges
end

function GraphShapeHist(block_approx::GraphHist{T,M}) where {T, M}
   max_shapes = get_num_blocks(block_approx) * (get_num_blocks(block_approx) + 1) ÷ 2
    return GraphShapeHist(max_shapes, block_approx)
end

function GraphShapeHist(number_shapes::Int,
        block_approx::GraphHist{T, M},
        X,
        k,
        n,
        max_num_shapes,
        indices_shape,
        indices_edges = CartesianIndex.((i, j) for j in 1:n for i in 1:n if i > j),
        algo = Yinyang()) where {
        T, M}
        Y = deepcopy(X)
    return GraphShapeHist!(Y, number_shapes, block_approx, X, k, n, max_num_shapes, indices_shape, indices_edges, algo)
end

function GraphShapeHist!(Y,number_shapes::Int,
        block_approx::GraphHist{T, M},
        X,
        k,
        n,
        max_num_shapes,
        indices_shape,
        indices_edges = CartesianIndex.((i, j) for j in 1:n for i in 1:n if i > j),
        algo = Yinyang()) where {
        T, M}
    if number_shapes >= max_num_shapes
        shape_labels = collect(1:max_num_shapes)
    else
        cluster_result = kmeans(algo, X, number_shapes)
        shape_labels = cluster_result.assignments
        for cluster in 1:number_shapes
            Y[:,shape_labels .== cluster] .= mean(X[:,shape_labels .== cluster], dims = 2)
        end
    end
    new_theta = Matrix{Array{T}}(undef, k, k)
    for (index,tuple_indices) in enumerate(indices_shape)
        new_theta[tuple_indices] = Y[:,index]
        new_theta[tuple_indices[2],tuple_indices[1]] = Y[:,index]
    end
    edge_labels = zeros(Int, n * (n - 1) ÷ 2)
    for i in 1:length(edge_labels)
        index_block = (block_approx.node_labels[indices_edges[i][1]], block_approx.node_labels[indices_edges[i][2]])
        index_block = index_block[1] ≥ index_block[2] ? index_block : (index_block[2], index_block[1])
        index_shape = findfirst(x -> Tuple(x) == index_block, indices_shape)
        edge_labels[i] = shape_labels[index_shape]
    end
    return GraphShapeHist(new_theta, block_approx.node_labels, edge_labels, block_approx.num_layers, number_shapes, shape_labels)
end


function GraphShapeHist(number_shapes::Int,block_approx::GraphHist{T,M}, algo = Yinyang()) where {T, M}
    X, k, n, max_num_shapes, indices_shape, indices_edges = make_x_matrix(block_approx)
    return GraphShapeHist(
        number_shapes, block_approx, X, k, n, max_num_shapes, indices_shape, indices_edges, algo)
end

function bic(g::GraphShapeHist{T,M}, A) where {T, M}
    n = length(g.node_labels)
    k = size(g.θ,1)
    indices = CartesianIndex.((i, j) for j in 1:k for i in 1:k if i ≥ j)
    number_parameters = (length(g.θ[1,1])-1)*length(unique(g.θ[indices]))
    return -2 * log_likelihood(g, A) + number_parameters * log(n * (n - 1) / 2)
end


function log_likelihood(g::Union{GraphShapeHist{T,M},GraphHist{T,M}}, A) where {T,M}
    assignment = Assignment(A, g.node_labels, GroupSize(g.node_labels))
    for i in 1:size(g.θ,1)
        for j in i:size(g.θ,2)
            assignment.estimated_theta[i,j,:] .= g.θ[i,j]
            assignment.estimated_theta[j,i,:] .= g.θ[i,j]
        end
    end
    return compute_log_likelihood(assignment)
end



function get_best_smoothed_estimator(g::GraphHist{T, M}, A; show_progress = true,
    n_min = get_num_blocks(g), n_max = binomial(get_num_blocks(g)+1, 2), steps = 1, max_iterations_stalled = 100,
        verbose = true, algo = Yinyang()) where {T, M}

    number_got_worse = 0
    if n_min > n_max
        n_min = n_max
    end

    best_bic = Inf
    best_g_shaped = nothing
    best_n_shape = n_max



    X, k, n, max_num_shapes, indices_shape, indices_edges = make_x_matrix(g)
    Y = deepcopy(X)
    max_shapes = min(n_max, length(unique(eachcol(X))))
    min_shapes = n_min
    if verbose
        println("Searching for the best number of shapes between $min_shapes and $max_shapes")
    end
    p = Progress(max_shapes; enabled = show_progress)
    n_shapes_iterator = min_shapes:steps:max_shapes
    #n_shapes_iterator = max_shapes:-steps:min_shapes
    bic_values = zeros(Union{Float64,Missing}, length(n_shapes_iterator))
    bic_values .= missing

    for (index, s) in enumerate(n_shapes_iterator)
        g_shaped = GraphShapeHist!(Y,
            s, g, X, k, n, max_num_shapes, indices_shape, indices_edges, algo)
        bic_value = bic(g_shaped, A)
        bic_values[index] = bic_value
        if bic_value < best_bic
            best_bic = bic_value
            best_g_shaped = deepcopy(g_shaped)
            best_n_shape = s
            number_got_worse = 0
        else
            number_got_worse += 1
        end
        if number_got_worse ≥ max_iterations_stalled
            break
        end
        next!(p; showvalues = [(:threshold, max_iterations_stalled),
                                (:stalled, number_got_worse),
                                (:best_n_shape, best_n_shape),
                                (:best_bic, best_bic),
                                (:current_n_shape, s),
                                (:current_bic, bic_value),
                                ]
            )
    end
    if show_progress
        finish!(p)
    end

    if verbose
        println("Best number of shapes: $best_n_shape")
    end

    return best_g_shaped, (n_shapes_iterator,bic_values)
end



# not type stable !
function get_moment_representation(g::GraphShapeHist{T, M}) where {T, M}
    if g.num_layers == 1
        return g.θ
    end
    moments = zeros(size(g.θ, 1), size(g.θ, 2), 2^g.num_layers - 1)
    transition = collect(kronecker([1 1; 0 1], g.num_layers))
    for i in 1:size(g.θ, 1)
        for j in 1:size(g.θ, 2)
            moments[i, j, :] .= (transition * g.θ[i, j])[2:end]
        end
    end
    indices_for_moments = [findall(x -> x == 1, _index_to_binary(e, g.num_layers))
                           for e in 2:size(g.θ[1, 1], 1)]
    return moments, indices_for_moments
end
