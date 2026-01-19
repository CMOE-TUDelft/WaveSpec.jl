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
using Statistics
using ..Integration
using RecipesBase

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
    low  = cdf(dist, Float64(a))
    high = cdf(dist, Float64(b))
    Z = high - low
    Z <= 0 && throw(ArgumentError("Zero probability mass in range [a, b]"))
    return TruncatedModel(dist, Float64(a), Float64(b), Z, low, high)
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
    try
        # Draw uniform u in [0, 1], scale to [CDF(a), CDF(b)]
        u = rand(rng) * tm.Z + tm.cdf_a
        return quantile(tm.dist, u)
    catch e
        # Case B: Fallback to Rejection Sampling
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

Statistics.mean(tm::TruncatedModel) = begin
    # Numerical integration for mean
    f(x) = x * pdf(tm, x)
    return IntegrateGaussQuad(f; a=tm.a, b=tm.b, order=4, n=100)
end

Statistics.var(tm::TruncatedModel) = begin
    μ = mean(tm)
    f(x) = (x - μ)^2 * pdf(tm, x)
    return IntegrateGaussQuad(f; a=tm.a, b=tm.b, order=4, n=100) - μ^2
end


# --- Plots and Visualisation ---
# Plot using light package RecipesBase:
#   1. import Plots in script
#   2. Call this function as plot(DiscreteSpectrum)

@recipe function f(s::TruncatedModel; n_points=500)
    # Plot Attributes
    title  := "Continuous Truncated Angular PDF"
    xlabel := "Angle θ (º)"
    ylabel := "Density S(θ) [m²s]"
    grid   := false

    # Get plot range 
    θ_range = range(s.a, s.b, length=n_points)

    @series begin
        label := false
        seriestype := :path
        fillrange := 0
        fillalpha := 0.2
        linecolor := :black
        θ_range, [pdf(s, θ) for θ in θ_range]
    end

end

end # module Truncation