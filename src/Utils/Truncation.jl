module Truncation

"""
Module Truncation
========================
This module defines a Truncated Distribution (TD) type by containing the necessary methods to truncate any given univariate distribution.
That is, the structure TruncatedModel takes a base distribution and truncates it between bounds [a, b], redefining the intrinsic univariate distribution methods:
    - probability density function (PDF)
    - cumulative distribution function (CDF)
    - random sampling
The resulting TD is a univariate distribution that is truncated to a specified interval [a, b], ensuring that all:
    - samples lie within this range
    - the distribution is properly normalized
"""

using Distributions
using Random

export TruncatedModel

# Define new concrete type: Truncated Distribution (TD)
struct TruncatedModel{D<:UnivariateDistribution, T<:Real} <: ContinuousUnivariateDistribution
    dist::D     # Base distribution
    a::T        # Lower bound
    b::T        # Upper bound
    Z::T        # Normalization constant Z = CDF(b) - CDF(a)
    cdf_a::T    # CDF(a)
    cdf_b::T    # CDF(b)
end

# Inner constructor (validates and precomputes constants)
function TruncatedModel(dist::UnivariateDistribution, a::Real, b::Real)
    low  = cdf(dist, a)
    high = cdf(dist, b)
    Z = high - low
    Z <= 0 && throw(ArgumentError("Zero probability mass in range [a, b]"))
    return TruncatedModel(dist, Float64(a), Float64(b), Float64(Z), Float64(low), Float64(high))
end

### Overload methods for truncated distribution: PDF, CDF and rand
# Truncated PDF
function Distributions.pdf(tm::TruncatedModel, x::Real)
    if x < tm.a || x > tm.b
        return 0.0
    else
        return pdf(tm.dist, x) / tm.Z
    end
end

# Truncated CDF
function Distributions.cdf(tm::TruncatedModel, x::Real)
    if x < tm.a
        return 0.0
    elseif x >= tm.b
        return 1.0
    else
        return (cdf(tm.dist, x) - tm.cdf_a) / tm.Z
    end
end

# Truncated Random Sampling
function Distributions.rand(rng::AbstractRNG, tm::TruncatedModel)

    # Case A: If the distribution supports quantile (Inverse Transform)
    if hasmethod(quantile, (typeof(tm.dist), Float64))
        # Draw uniform u in [0, 1], scale to [CDF(a), CDF(b)]
        u = rand(rng) * tm.Z + tm.cdf_a
        return quantile(tm.dist, u)

    # Case B: Fallback to Rejection Sampling
    else
        while true
            x = rand(rng, tm.dist)
            if tm.a <= x <= tm.b
                return x
            end
        end
    end
end

# 2. Vectorized method for multiple samples (e.g., rand(rng, tm, 10))
function Distributions.rand(rng::AbstractRNG, tm::TruncatedModel, n::Int)
    # Pre-allocate the vector for efficiency
    samples = Vector{Float64}(undef, n)
    for i in 1:n
        samples[i] = rand(rng, tm) # This calls your method above
    end
    return samples
end

# convenience rand when rng::AbstractRNG random generator type is not passed as an argument
Distributions.rand(tm::TruncatedModel) = rand(Random.GLOBAL_RNG, tm)

end # module Truncation