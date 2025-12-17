
"""
    DonelanBannerDistribution(μ, β)

The Donelan-Banner spreading distribution over (-∞, ∞), where:
    - pdf: D(θ) = (β/2) * sech²(β(θ - μ))
    - cdf: F(θ) = 0.5 * (1 + tanh(β(θ - μ)))
    - quantile: F⁻¹(u) = μ + atanh(2u - 1) / β
    - `μ`: Mean direction (radians).
    - `β`: Spreading parameter (β > 0). 
"""

export DonelanBannerDistribution

struct DonelanBannerDistribution{T<:Real} <: ContinuousUnivariateDistribution
    μ::T               # Mean direction
    β::T               # Spreading parameter
end

# Inner constructor with validation
function DonelanBannerDistribution(μ::T, β::T) where {T<:Real}
    β > 0 || throw(ArgumentError("Spreading parameter β must be positive"))
    return DonelanBannerDistribution{T}(μ, β)
end

# Outer constructor
DonelanBannerDistribution(μ::Real, β::Real) = DonelanBannerDistribution(promote(μ, β)...)

# --- PDF ---
# D(θ) = (β/2) * sech²(β(θ - μ))
function Distributions.pdf(d::DonelanBannerDistribution, x::Real)
    return (d.β / 2.0) * (sech(d.β * (x - d.μ))^2)
end

# --- CDF ---
# F(x) = 0.5 * (1 + tanh(β(x - μ)))
function Distributions.cdf(d::DonelanBannerDistribution, x::Real)
    return 0.5 * (1.0 + tanh(d.β * (x - d.μ)))
end

# --- Quantile (Inverse CDF) ---
# This allows the Truncation module to perform Inverse Transform Sampling
# quantile: F⁻¹(u) = μ + atanh(2u - 1) / β
function Distributions.quantile(d::DonelanBannerDistribution, u::Real)
    # Check bounds for safety
    (u <= 0) && return -Inf
    (u >= 1) && return Inf
    return d.μ + atanh(2.0 * u - 1.0) / d.β
end

# --- Sampling ---
function Distributions.rand(rng::AbstractRNG, d::DonelanBannerDistribution)
    return quantile(d, rand(rng))
end

# convenience rand when rng::AbstractRNG random generator type is not passed as an argument
Distributions.rand(d::DonelanBannerDistribution) = rand(Random.GLOBAL_RNG, d)
