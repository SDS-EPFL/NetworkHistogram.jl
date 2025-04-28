# Inspired by Discretizer.jl but with the fast decoding function and built-in
# convention for discretizing continuous distributions.
abstract type Discretizer end

function encode(d::Discretizer, x::AbstractArray{<:Real})
    return [encode(d, u) for u in x]
end

function decode(d::Discretizer, x::AbstractArray{<:Real})
    return [decode(d, u) for u in x]
end

"""
Uniformly discretizes a continuous distribution into a fixed number of bins of equal width.
"""
struct RegularDiscretizer{F, T, L} <: Discretizer
    n_bins::Int
    lower_bound::F
    upper_bound::F
    bin_labels::MVector{L, T}
    bin_width::F
end

function RegularDiscretizer(
        n_bins::Int, lower_bound::F, upper_bound::F) where {F}
    if !isfinite(lower_bound) || !isfinite(upper_bound)
        throw(ArgumentError("RegularDiscretizer requires finite lower and upper bounds."))
    end
    bin_width = (upper_bound - lower_bound) / n_bins
    return RegularDiscretizer(
        n_bins, lower_bound, upper_bound, MVector{n_bins}(1:n_bins), bin_width
    )
end

function support_encoding(d::RegularDiscretizer, x::Real)
    return d.lower_bound <= x <= d.upper_bound
end

function minimum(d::RegularDiscretizer)
    return d.lower_bound
end

function maximum(d::RegularDiscretizer)
    return d.upper_bound
end

function encode(d::RegularDiscretizer, x::Real)
    if x >= d.upper_bound
        return d.n_bins
    end
    return d.bin_labels[convert(Int, div(x - d.lower_bound, d.bin_width) + 1)]
end

function _decode_randomly(
        rng::Random.AbstractRNG, d::RegularDiscretizer, bin::Int)
    hi, lo = decode(d, bin)
    return lo + (hi - lo) * rand(rng)
end

function binwidth(d::RegularDiscretizer)
    return d.bin_width
end

function decode(d::RegularDiscretizer, bin::Int)
    return (d.lower_bound + (bin - 1) * d.bin_width,
        d.lower_bound + bin * d.bin_width)
end

function encode(d::RegularDiscretizer, x::AbstractArray{Real})
    return [encode(d, u) for u in x]
end

function decode(d::RegularDiscretizer, x::AbstractArray{Real})
    return [decode(d, u) for u in x]
end

function nlabels(d::RegularDiscretizer)
    return d.n_bins
end

non_zero_labels_counts(d::RegularDiscretizer) = nlabels(d)

"""
Maps a set of categories to a set of bins
"""
struct CategoryDiscretizer{F, T}
    cat_to_bin::Dict{F, T}
    bin_to_cat::Dict{T, F}
    min_label::T
    max_label::T
end

function CategoryDiscretizer(cat_to_bin::Dict, bin_to_cat::Dict)
    min_label = minimum(keys(bin_to_cat))
    max_label = maximum(keys(bin_to_cat))
    return CategoryDiscretizer(cat_to_bin, bin_to_cat, min_label, max_label)
end

function support_encoding(d::CategoryDiscretizer, x)
    return haskey(d.cat_to_bin, x)
end

function encode(d::CategoryDiscretizer, x)
    return d.cat_to_bin[x]
end

function decode(d::CategoryDiscretizer, label)
    return d.bin_to_cat[label]
end

function nlabels(d::CategoryDiscretizer)
    return length(d.bin_to_cat)
end

function binwidth(d::CategoryDiscretizer{F, T}, x::T) where {F, T}
    return length(d.bin_to_cat[x])
end

function non_zero_labels_counts(d::CategoryDiscretizer)
    if 0 ∈ keys(d.bin_to_cat)
        return length(d.bin_to_cat) - 1
    else
        return length(d.bin_to_cat)
    end
end

function minimum(d::CategoryDiscretizer)
    return d.min_label
end

function maximum(d::CategoryDiscretizer)
    return d.max_label
end

"""
Uniformly discretizes a continuous distribution into a fixed number of bins of equal width,
with additional bins for missing or special values.
"""
struct HybridDiscretizer{F, T, L} <: Discretizer
    lin::RegularDiscretizer{F, T, L}
    cat::CategoryDiscretizer{F, T}
end

# change so that atoms can be packed together if wanted
function HybridDiscretizer(n_bins, lower_bound, upper_bound, atoms)
    cat_to_bin = Dict(a => n_bins + i for (i, a) in enumerate(atoms))
    bin_to_cat = Dict(n_bins + i => a for (i, a) in enumerate(atoms))
    bin_width = (upper_bound - lower_bound) / n_bins
    return HybridDiscretizer(
        RegularDiscretizer{typeof(bin_width), Int, n_bins}(
            n_bins, lower_bound, upper_bound, MVector{n_bins}(1:n_bins),
            (upper_bound - lower_bound) / n_bins),
        CategoryDiscretizer(cat_to_bin, bin_to_cat)
    )
end

function DiscretizerZeroToZero(n_bins, lower_bound, upper_bound)
    cat_to_bin = Dict([0.0 => 0])
    bin_to_cat = Dict([0 => 0.0])
    bin_width = (upper_bound - lower_bound) / n_bins
    return HybridDiscretizer(
        RegularDiscretizer{typeof(bin_width), Int, n_bins}(
            n_bins, lower_bound, upper_bound, MVector{n_bins}(1:n_bins),
            (upper_bound - lower_bound) / n_bins),
        CategoryDiscretizer(cat_to_bin, bin_to_cat)
    )
end

function support_encoding(d::HybridDiscretizer, x)
    return support_encoding(d.lin, x) || support_encoding(d.cat, x)
end

function minimum(d::HybridDiscretizer)
    return min(minimum(d.lin), minimum(d.cat))
end

function maximum(d::HybridDiscretizer)
    return max(maximum(d.lin), maximum(d.cat))
end

function nlabels(d::HybridDiscretizer)
    return nlabels(d.lin) + nlabels(d.cat)
end

function non_zero_labels_counts(d::HybridDiscretizer)
    return non_zero_labels_counts(d.lin) + non_zero_labels_counts(d.cat)
end

binwidth(d::HybridDiscretizer) = binwidth(d.lin)

function binwidth(d::HybridDiscretizer, bin)
    if haskey(d.cat.cat_to_bin, bin)
        return binwidth(d.cat, bin)
    else
        return binwidth(d.lin)
    end
end

function encode(d::HybridDiscretizer, x::Real)
    if haskey(d.cat.cat_to_bin, x)
        return encode(d.cat, x)
    else
        return encode(d.lin, x)
    end
end

function decode(d::HybridDiscretizer, bin::Int)
    if haskey(d.cat.bin_to_cat, bin)
        return decode(d.cat, bin)
    else
        return decode(d.lin, bin)
    end
end

function _decode_randomly(
        rng::Random.AbstractRNG, d::HybridDiscretizer, bin::Int)
    if haskey(d.cat.bin_to_cat, bin)
        return decode(d.cat, bin)
    else
        return _decode_randomly(rng, d.lin, bin)
    end
end

function auto_nbins(data)
    binwidth = 2iqr(data) / cbrt(n)
    lo, hi = extrema(data)
    nbins_fd = ceil(Int, (hi - lo) / binwidth)
    nbins_sturges = ceil(Int, log(2, n)) + 1
    nbins = max(nbins_fd, nbins_sturges)
    return nbins
end

function progress_in_bin(d::CategoryDiscretizer, x::Real, bin)
    return one(x)
end

function progress_in_bin(d::RegularDiscretizer, x::Real, bin)
    lo, hi = decode(d, bin)
    return (x - lo) / (hi - lo)
end

function progress_in_bin(d::HybridDiscretizer, x::Real, bin)
    if haskey(d.cat.bin_to_cat, bin)
        return progress_in_bin(d.cat, x, bin)
    else
        return progress_in_bin(d.lin, x, bin)
    end
end
