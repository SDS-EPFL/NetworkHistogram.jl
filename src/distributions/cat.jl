const Cat{M, T} = Categorical{T, SVector{M, T}}

function Cat(p::SVector{M, T}) where {M, T}
    return Categorical(p)
end

function Base.show(io::IO, c::Cat)
    Base.print(io, "Cat($(c.p))")
end

num_categories(::Type{Cat{M, T}}) where {M, T} = M
num_categories(::Cat{M, T}) where {M, T} = M
zero(c::Cat{M, T}) where {M, T} = Cat(ones(typeof(c.p)) ./ M)
sample(c::Cat) = rand(c)
function fit(c::Cat{M, T}, xs::AbstractVector{Int}) where {M, T}
    total = length(xs)
    if total == 0
        return zero(c)
    end
    return Cat(SVector{M}(counts(xs, M) ./ total))
end

function fit(::Cat{M, T}, x::Int) where {M, T}
    ps = zeros(T, M)
    ps[x] = one(T)
    return Cat(SVector{M}(ps))
end

function _xlogy(x, y)
    if x == 0
        return zero(y)
    end
    return x * log(y)
end

function logpdf_cat(p::AbstractVector, obs::Int)
    return log(p[obs])
end

# Efficient log-likelihood computation for categorical observations
# Uses xlogy(x,y) = x*log(y) which handles edge cases properly
function logpdf_cat(p::AbstractVector, count_observed::AbstractVector)
    return sum(_xlogy.(count_observed, p))
end

distance(c1::Cat{M, V}, c2::Cat{M, V}) where {M, V} = sum(abs.(c1.p .- c2.p))

function get_ref_dist(dist::Categorical, ::Val{true})
    return Dist(Cat(SVector{ncategories(dist) + 1}(0.0, dist.p...)))
end

function get_ref_dist(dist::Categorical, ::Val{false})
    return Dist(Cat(SVector{ncategories(dist)}(dist.p)))
end

_fast_compressed_obs(d::Categorical, x::Int, ::Val{true}) = x + one(x)
_fast_compressed_obs(d::Categorical, x::Int, ::Val{false}) = x

function tv_distance(c1::Cat, c2::Cat)
    return sum(abs.(c1.p .- c2.p)) / 2
end

function l2_distance(c1::Cat, c2::Cat)
    return sqrt(sum((c1.p .- c2.p) .^ 2))
end
