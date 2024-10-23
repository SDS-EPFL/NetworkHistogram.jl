mutable struct CategoricalData{M, F, C}
    counts::Matrix{Int}
    realized::Matrix{MVector{M, Int}}
    estimated_theta::Matrix{MVector{M, F}}
    A::Matrix{C} # possible use of CategoricalArrays.jl ?
    log_likelihood::F
end

const CategoricalAssignment{T, M, F, C} = Assignment{
    T, CategoricalData{M, F, C}}
const CategoricalInitRule{S, F} = InitRule{S, Val{CategoricalData}}

function CategoricalAssignment(
        g, group_size::GroupSize, node_labels::Vector{Int})
    categorical_data = make_categorical_data(g, node_labels, group_size)
    return Assignment(group_size, node_labels, categorical_data)
end

function make_assignment(g, h, init_rule::CategoricalInitRule)
    group_size,
    node_labels = initialize_node_labels(
        g, h, init_rule.starting_assignment_rule)
    a = CategoricalAssignment(deepcopy(g), group_size, node_labels)
    @show a.additional_data.log_likelihood
    ll_test = force_recompute_ll(a, g)
    @show ll_test
    return a
end

function make_categorical_data(g, node_labels, group_size)
    number_groups = length(group_size)
    A, num_categories = categorical_matrix(g)
    counts = zeros(Int, number_groups, number_groups)
    realized = [MVector{num_categories}(zeros(Int, num_categories))
                for _ in 1:number_groups, _ in 1:number_groups]

    _count_cat_occurences!(counts, realized, g, Assignment(group_size, node_labels))

    estimated_theta = realized ./ counts
    ll = compute_log_likelihood(estimated_theta, realized)
    return CategoricalData(counts, realized, estimated_theta, A, ll)
end

function _count_cat_occurences!(counts, realized, g, a_dummy)
    @inbounds for k in 1:number_groups(a_dummy)
        for l in k:number_groups(a_dummy)
            counts_dict = StatsBase.countmap(get_obs.(
                Ref(g), get_edge_indices(a_dummy, k, l)))
            total = 0
            for (m, v) in counts_dict
                realized[k, l][m] = v
                realized[l, k][m] = v
                total += v
            end
            counts[k, l] = total
            counts[l, k] = total
        end
    end
end


function recount_occurences!(a)
    _count_cat_occurences!(a.additional_data.counts, a.additional_data.realized, a.additional_data.A, a)
    return nothing
end

function compute_log_likelihood(
        estimated_theta::AbstractMatrix{MVector{M, T}}, realized::AbstractMatrix{F}) where {
        M, T, F}
    loglik = zero(T)
    number_groups = size(estimated_theta, 1)
    @inbounds for j in 1:number_groups
        for i in j:number_groups
            for m in 1:length(estimated_theta[i, j])
                if realized[i, j][m] != 0
                    loglik += realized[i, j][m] * log(estimated_theta[i, j][m])
                end
            end
            #loglik += sum(log.(estimated_theta[i, j]) .* realized[i, j])
            #loglik += sum(xlogy.(realized[i,j], estimated_theta[i, j]) )
        end
    end
    return loglik
end

function categorical_matrix(A::CategoricalMatrix)
    @info "Converting CategoricalMatrix to matrix"
    categories = levels(A)
    return levelcode.(recode(
        A, [l => i for (i, l) in enumerate(categories)]..., missing => 0))
end

# to update, just for test now
function categorical_matrix(A::AbstractMatrix{Int})
    min_A = minimum(A)
    if min_A > 1
        A_inter = A .- min_A .+ 1
    else
        A_inter = copy(A)
    end
    for i in 1:size(A_inter, 1)
        A_inter[i, i] = 0
    end
    return A_inter
end

function categorical_matrix(g::Observations)
    return categorical_matrix(g.graph), length(support(g.dist_ref))
end

function loglikelihood(a::CategoricalAssignment, g::Observations)
    return a.additional_data.log_likelihood
end

function force_recompute_ll(a::CategoricalAssignment, g::Observations)
    a_simple = Assignment(a.group_size, a.node_labels)
    return loglikelihood(a_simple, g)
end

include("swap.jl")
