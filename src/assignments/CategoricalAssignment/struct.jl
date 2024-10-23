mutable struct CategoricalData{F, C}
    counts::Matrix{Int}
    realized::Array{Int, 3}
    estimated_theta::Array{F, 3}
    A::Matrix{C} # possible use of CategoricalArrays.jl ?
    log_likelihood::F # need to remove this type
end

const CategoricalAssignment{T, F, C} = Assignment{
    T, CategoricalData{F, C}}
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
    a = CategoricalAssignment(g, group_size, node_labels)
    return a
end

function make_categorical_data(g, node_labels, group_size)
    number_groups = length(group_size)
    A, num_categories = categorical_matrix(g)
    counts = zeros(Int, number_groups, number_groups)
    realized = zeros(Int, num_categories, number_groups, number_groups)
    estimated_theta = zeros(
        Float64, num_categories, number_groups, number_groups)

    _count_cat_occurences!(
        counts, realized, g, Assignment(group_size, node_labels))

    _fast_div!(estimated_theta, realized, counts)

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
                realized[m, k, l] = v
                realized[m, l, k] = v
                total += v
            end
            counts[k, l] = total
            counts[l, k] = total
        end
    end
end

function recount_occurences!(a)
    _count_cat_occurences!(
        a.additional_data.counts, a.additional_data.realized, a.additional_data.A, a)
    return nothing
end

function compute_log_likelihood(
        estimated_theta::Array{T, 3}, realized::Array{F, 3}) where {
        T, F}
    loglik = zero(T)
    number_groups = size(estimated_theta, 2)
    number_decorations = size(estimated_theta, 1)
    @inbounds for j in 1:number_groups
        for i in j:number_groups
            for m in 1:number_decorations
                if realized[m, i, j] != 0
                    loglik += realized[m, i, j] * log(estimated_theta[m, i, j])
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
