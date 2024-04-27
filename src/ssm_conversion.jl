struct GraphShapeHist{T,M}
    θ::Array{T, M}
    node_labels::Vector{Int}
    edge_labels::Vector{Int}
    num_layers::Int
end

# only work for 1D ....
function GraphShapeHist(number_shapes::Int,block_approx::GraphHist{T,M}) where {T, M}
    k = get_num_blocks(block_approx)
    n = length(block_approx.node_labels)
    num_shapes = k*(k+1)÷2
    X = Matrix{T}(undef, get_num_params_per_block(block_approx), num_shapes)
    mapping_block_to_shape = Dict{Tuple{Int,Int},Int}()
    mapping_shape_to_block = Dict{Int,Array{Tuple{Int,Int}}}()
    index = 1
    for i in 1:k
        for j in i:k
            X[:,index] = block_approx.θ[i,j,:]
            mapping_block_to_shape[(i,j)] = index
            if haskey(mapping_shape_to_block,index)
                push!(mapping_shape_to_block[index],(i,j))
            else
                mapping_shape_to_block[index] = [(i,j)]
            end
            index += 1
        end
    end
    cluster_result = kmeans(X, number_shapes; maxiter = 1000)
    shape_labels = cluster_result.assignments
    shape_values = Matrix{T}(undef, number_shapes, get_num_params_per_block(block_approx))
    for i in 1:nclusters(cluster_result)
        shape_values[i,:] .= mean(X[:,shape_labels .== i], dims = 2)
    end
    new_theta = Matrix{Array{T}}(undef, k, k)
    edge_labels = zeros(Int, n*(n-1)÷2)
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
    return GraphShapeHist(new_theta, block_approx.node_labels, edge_labels, block_approx.num_layers)
end
