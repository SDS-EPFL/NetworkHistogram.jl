struct GraphShapeHist{T,M}
    θ::Array{T, M}
    node_labels::Vector{Int}
    edge_labels::Vector{Int}
    num_layers::Int
    num_shapes::Int
end

function GraphShapeHist(number_shapes::Int,block_approx::GraphHist{T,M}) where {T, M}
    k = get_num_blocks(block_approx)
    n = length(block_approx.node_labels)
    num_shapes = k*(k+1)÷2
    X = Matrix{T}(undef, get_num_params_per_block(block_approx), num_shapes)
    mapping_block_to_shape = Dict{Tuple{Int,Int},Int}()
    mapping_shape_to_block = Dict{Int,Array{Tuple{Int,Int}}}()
    index = 1
    edge_labels = zeros(Int, n * (n - 1) ÷ 2)
    new_theta = Matrix{Array{Float64, 1}}(undef, k, k)
    for i in 1:k
        for j in i:k
            inter = block_approx.θ[i,j,:]
            new_theta[i,j] = inter
            X[:,index] = inter
            mapping_block_to_shape[(i,j)] = index
            if haskey(mapping_shape_to_block,index)
                push!(mapping_shape_to_block[index],(i,j))
            else
                mapping_shape_to_block[index] = [(i,j)]
            end
            index += 1
        end
    end
    if number_shapes == num_shapes # no need to cluster
        for i in 1:n
            for j in i+1:n
                label_i = block_approx.node_labels[i]
                label_j = block_approx.node_labels[j]
                if label_i > label_j
                    label_i,label_j = label_j,label_i
                end
                index = mapping_block_to_shape[(label_i, label_j)]
                edge_labels[index] = mapping_block_to_shape[(label_i, label_j)]
            end
        end
        return GraphShapeHist(new_theta, block_approx.node_labels,
            edge_labels, block_approx.num_layers, number_shapes)
    end

    cluster_result = kmeans(X, number_shapes; maxiter = 1000)
    shape_labels = cluster_result.assignments
    shape_values = Matrix{T}(undef, number_shapes, get_num_params_per_block(block_approx))
    for i in 1:nclusters(cluster_result)
        shape_values[i,:] .= mean(X[:,shape_labels .== i], dims = 2)
    end

    for i in 1:k
        for j in i:k
            shape_index = shape_labels[mapping_block_to_shape[(i,j)]]
            new_theta[i,j] = shape_values[shape_index,:]
            new_theta[j,i] = shape_values[shape_index,:]
        end
    end

    index = 1
    for i in 1:n
        for j in i+1:n
            label_i = block_approx.node_labels[i]
            label_j = block_approx.node_labels[j]
            if label_i > label_j
                label_i,label_j = label_j,label_i
            end
            edge_labels[index] = shape_labels[mapping_block_to_shape[(label_i,label_j)]]
            index += 1
        end
    end
    if size(new_theta[1,1],1) == 1
        inter_theta = zeros(T, k, k)
        for i in 1:k
            for j in i:k
                inter_theta[i,j] = new_theta[i,j][1]
                inter_theta[j,i] = new_theta[i,j][1]
            end
        end
        new_theta = inter_theta
    end
    return GraphShapeHist(new_theta, block_approx.node_labels, edge_labels, block_approx.num_layers, number_shapes)
end

function bic(g::GraphShapeHist{T,M}, A) where {T, M}
    n = length(g.node_labels)
    number_parameters = length(unique(g.edge_labels))
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



function get_best_smoothed_estimator(g::GraphHist{T,M}, A) where {T,M}
    group_number = length(unique(g.node_labels))
    max_shapes = group_number*(group_number+1)÷2
    best_bic = Inf
    best_g_shaped = nothing
    bic_values = zeros(max_shapes)
    @showprogress for s in max_shapes:-1:1
        g_shaped = GraphShapeHist(s, g)
        bic_value = bic(g_shaped, A)
        bic_values[s] = bic_value
        @debug "BIC for $s shapes: ", bic_value
        if bic_value < best_bic
            best_bic = bic_value
            best_g_shaped = g_shaped
        end
    end
    return best_g_shaped, bic_values
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
