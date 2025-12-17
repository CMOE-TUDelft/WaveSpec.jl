
"""
    CosinePowerDistribution(μ, n)

This module contains the Cosine Power distribution defined over [μ - π, μ + π] with:
    - `μ`: Mean direction (radians)
    - `n`: Power exponent (n > 0). 

Special cases:
    - n = 0: Uniform distribution
    - n = 2: Cosine Squared distribution
"""

export CosinePowerDistribution

struct CosinePowerDistribution{T<:Real} <: ContinuousUnivariateDistribution
    μ::T               # Mean direction
    n::T               # Power exponent
    norm_factor::T     # Normalization constant
    
    # Inner constructor with validation
    function CosinePowerDistribution(μ::T, n::T) where {T<:Real}
        # Validate parameters
        n > 0 || throw(ArgumentError("Power n must be positive"))

        # The integral of cos(x/2)^n from -π to π is 2^(n+1) * Beta(n/2 + 1, n/2 + 1)
        # However, for the standard [-π, π] range, the factor is simply:
        # B(1/2, n/2 + 1)
        norm_factor = 1.0 / (2 * beta(0.5, n/2 + 1))

        return new{T}(μ, n, norm_factor)      
    end
end

# Shorthand constructor for different types
CosinePowerDistribution(μ::Real, n::Real) = CosinePowerDistribution(promote(μ, n)...)

# --- PDF ---
function Distributions.pdf(d::CosinePowerDistribution, θ::Real)
    # Map input angle back to [-π, π]
    diff = mod(θ - d.μ + π, 2π) - π
    
    # cos(diff/2) is always non-negative in [-π, π]
    return d.norm_factor * cos(diff/2)^d.n
end

# --- CDF ---
function Distributions.cdf(d::CosinePowerDistribution, θ::Real)
    # Map input angle back to [-π, π]
    diff = mod(θ - d.μ + π, 2π) - π  
    
    # The CDF of cos(x/2)^n relates to the Incomplete Beta Function
    # Let x = sin(diff/2)^2. Then we use the Regularized Incomplete Beta.
    # This is much faster than numerical integration.
    
    val = 0.5 * beta_inc(0.5, d.n/2 + 1, sin(diff/2)^2)[1]
    
    return (diff >= 0) ? (0.5 + val) : (0.5 - val)
end

# --- Sampling ---
# Since we have the CDF (via Beta), but not an easy Inverse CDF, 
# Rejection Sampling is the most reliable and easiest to implement.
function Distributions.rand(rng::AbstractRNG, d::CosinePowerDistribution)
    # Peak value of PDF is at θ = μ, where cos(0)=1
    max_pdf = 1.0 / (2 * beta(0.5, d.n/2 + 1))
    
    while true
        θ = (rand(rng) * 2π) - π + d.μ
        y = rand(rng) * max_pdf
        if y <= pdf(d, θ)
            return θ
        end
    end
end

# convenience rand when rng::AbstractRNG random generator type is not passed as an argument
Distributions.rand(d::CosinePowerDistribution) = rand(Random.GLOBAL_RNG, d)
