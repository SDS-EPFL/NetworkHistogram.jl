"""
    Dist{D}

A wrapper for distributions that tracks aggregation statistics.

This type wraps a distribution `D` and maintains a count of how many observations
have been aggregated into it. This is essential for the network histogram algorithm
which needs to efficiently update distributions as nodes move between groups.

# Fields
- `dist::D`: The underlying distribution
- `counts::Int`: Number of observations aggregated into this distribution (must be ≥ 0)

# Type Parameters
- `D`: The type of the underlying distribution (e.g., Bernoulli, Categorical, etc.)

# Constructors
```julia
# With explicit count
Dist(distribution, counts::Int)

# Single observation (count = 1)
Dist(distribution)
```

# Examples
```julia
# Wrap a Bernoulli distribution
d = Dist(Bernoulli(0.5))

# Create a zero distribution
d_zero = zero(d)

# Add observations
d_updated = add_to(d, Bernoulli(0.7))

# Remove observations
d_reduced = remove_from(d_updated, Bernoulli(0.7))
```

See also: [`add_to`](@ref), [`remove_from`](@ref), [`zero`](@ref)
"""
struct Dist{D}
    dist::D
    counts::Int
    function Dist(d, counts::Int)
        if counts < 0
            throw(ArgumentError("Counts ($counts) cannot be negative"))
        end
        new{typeof(d)}(d, counts)
    end
end

function Base.show(io::IO, d::Dist)
    print(io, "$(d.dist)")
end

"""
    Dist(d)

Create a Dist with a single observation (count = 1).
"""
Dist(d) = Dist(d, 1)

"""
    zero(d::Dist)

Create a zero-initialized distribution with 0 counts.
"""
zero(d::Dist) = Dist(zero(d.dist), 0)

Base.broadcastable(x::Dist) = Ref(x)

"""
    add_to(avgdist::Dist{D}, dist::D) where {D}

Add a new observation to an aggregated distribution.

Updates the distribution parameters using weighted averaging based on the count.
The new observation has weight 1/(counts+1) and the existing distribution has
weight counts/(counts+1).

# Arguments
- `avgdist::Dist{D}`: The current aggregated distribution
- `dist::D`: The new distribution to add

# Returns
- `Dist{D}`: Updated distribution with incremented count

# Example
```julia
d = Dist(Bernoulli(0.5), 2)  # 2 observations with mean 0.5
d_new = add_to(d, Bernoulli(0.8))  # Add observation with value 0.8
# Result: Dist with 3 observations and mean (2*0.5 + 1*0.8)/3 ≈ 0.6
```
"""
function add_to(avgdist::Dist{D}, dist::D) where {D}
    inner_dist = agg_params(
        avgdist.dist, dist, avgdist.counts / (avgdist.counts + 1),
        1 / (avgdist.counts + 1))
    return Dist(inner_dist, avgdist.counts + 1)
end

"""
    remove_from(avgdist::Dist{D}, dist::D) where {D}

Remove an observation from an aggregated distribution.

Updates the distribution parameters by removing the contribution of `dist` from
the aggregate, using appropriate weight adjustments.

# Arguments
- `avgdist::Dist{D}`: The current aggregated distribution
- `dist::D`: The distribution to remove

# Returns
- `Dist{D}`: Updated distribution with decremented count

# Note
Throws an error if attempting to remove from a distribution with 0 counts.
"""
function remove_from(avgdist::Dist{D}, dist::D) where {D}
    if avgdist.counts <= 0
        error("Cannot remove from a distribution with 0 counts")
    end
    return Dist(
        agg_params(
            avgdist.dist, dist, avgdist.counts / max(1, (avgdist.counts - 1)),
            -1 / max(1, (avgdist.counts - 1))),
        avgdist.counts - 1)
end

"""
    add_to(avgdist::Dist{D}, dist::Dist{D}) where {D}

Add two Dist objects together, properly accounting for their counts.

# Arguments
- `avgdist::Dist{D}`: First distribution
- `dist::Dist{D}`: Second distribution to add

# Returns
- `Dist{D}`: Combined distribution with summed counts
"""
function add_to(avgdist::Dist{D}, dist::Dist{D}) where {D}
    Dist(
        agg_params(
            avgdist.dist, dist.dist, avgdist.counts /
                                     (avgdist.counts + dist.counts),
            dist.counts / (avgdist.counts + dist.counts)),
        avgdist.counts + dist.counts)
end

"""
    remove_from(avgdist::Dist, dist::Dist)

Remove one Dist from another, properly accounting for their counts.
"""
function remove_from(avgdist::Dist, dist::Dist)
    Dist(
        agg_params(
            avgdist.dist, dist.dist,
            avgdist.counts / max(1, (avgdist.counts - dist.counts)),
            -dist.counts / max(1, (avgdist.counts - dist.counts))),
        avgdist.counts - dist.counts)
end

"""
    _fast_compressed_obs(d, x, zero_inflated)

Compress observations for efficient storage and computation.

By default, returns `x` unchanged. Distributions can override this to implement
custom compression strategies.
"""
_fast_compressed_obs(d, x, zero_inflated) = x

# Delegate common operations to the underlying distribution
for f in [:logpdf, :sample, :distance, :eltype, :params, :_fast_compressed_obs]
    @eval $f(d::Dist, args...) = $f(d.dist, args...)
end

"""
    fit(d::Dist, x)

Fit the underlying distribution to observation(s) `x`, preserving the count.
"""
fit(d::Dist, x) = Dist(fit(d.dist, x), d.counts)

"""
    loglikelihood(d::Dist, x)

Compute the log-likelihood of observation(s) `x` under distribution `d`.

# Returns
- `Float64`: Sum of log-probabilities, or 0.0 if x is empty
"""
loglikelihood(d::Dist, x) = isempty(x) ? 0.0 : sum(logpdf(d, y) for y in x)

"""
    unwrap(d::Dist)

Extract the underlying distribution from a Dist wrapper.
"""
unwrap(d::Dist) = d.dist

Base.promote_rule(::Type{Dist{D}}, ::Type{D}) where {D} = D
Base.convert(::Type{D}, d::Dist{D}) where {D} = d.dist

"""
    Bernoulli{T <: Real}

A simple Bernoulli distribution for binary (0/1) edges.

# Fields
- `p::T`: Success probability (probability of edge = 1)

# Example
```julia
b = Bernoulli(0.3)  # 30% chance of edge
edge = sample(b)     # Returns true or false
ll = logpdf(b, true) # Log probability of observing an edge
```

# Interface Requirements
For a distribution to work with NetworkHistogram, it must implement:
- `zero(d)`: Return a zero-initialized distribution
- `agg_params(d1, d2, w1, w2)`: Aggregate two distributions with weights
- `fit(d, x)`: Fit distribution to observation(s)
- `distance(d1, d2)`: Distance metric between distributions
- `logpdf(d, x)`: Log probability density/mass function
- `params(d)`: Return tuple of parameters
- `eltype(d)`: Return element type
- `sample(d)`: Generate a random sample
"""
struct Bernoulli{T <: Real}
    p::T
    function Bernoulli(p::T) where {T <: Real}
        if isnan(p) || isinf(p)
            throw(ArgumentError("Bernoulli parameter p=$p must be finite"))
        end
        if !(0 <= p <= 1)
            throw(ArgumentError("Bernoulli parameter p=$p must be in [0, 1]"))
        end
        new{T}(p)
    end
end

zero(d::Bernoulli) = Bernoulli(zero(d.p))
zero(::Type{Bernoulli{T}}) where {T} = Bernoulli(zero(T))
function agg_params(d1::Bernoulli, d2::Bernoulli, w1, w2)
    p = w1 * d1.p + w2 * d2.p
    # Clamp to [0, 1] to handle floating-point arithmetic errors
    p = clamp(p, 0.0, 1.0)
    Bernoulli(p)
end
fit(::Bernoulli, x) = Bernoulli(mean(x))
distance(d1::Bernoulli, d2::Bernoulli) = abs(d1.p - d2.p)
logpdf(d::Bernoulli, x) = log(d.p * x + (1 - d.p) * (1 - x))
params(d::Bernoulli) = (d.p,)
eltype(d::Bernoulli) = Bool
sample(d::Bernoulli) = rand() <= d.p
