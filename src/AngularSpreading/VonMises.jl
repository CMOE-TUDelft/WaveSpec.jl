
"""
    VonMisesDistribution(μ, κ)

The Von Mises (Circular Normal) distribution.
PDF: D(θ) = exp(κ * cos(θ - μ)) / (2π * I₀(κ))
- `μ`: Mean direction (radians).
- `κ`: Concentration parameter (κ > 0). κ=0 is uniform; κ→∞ is Normal.
"""

export VonMisesDistribution

struct VonMisesDistribution{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    κ::T
    norm_factor::T # Pre-compute log of the normalization factor for stability
end

function VonMisesDistribution(μ::T, κ::T) where {T<:Real}
    κ >= 0 || throw(ArgumentError("Concentration κ must be non-negative"))
    
    # Normalization factor C = 2π * I₀(κ)
    # Using besselix(0, κ) = exp(-|κ|) * I₀(κ) to prevent overflow at high κ
    C = 1.0 / (2π * besselix(0, κ) * exp(κ))

    return VonMisesDistribution{T}(μ, κ, C)
end

VonMisesDistribution(μ::Real, κ::Real) = VonMisesDistribution(promote(μ, κ)...)

# --- PDF ---
function Distributions.pdf(d::VonMisesDistribution, x::Real)
    return d.norm_factor * exp(d.κ * cos(x - d.μ)) 
end

# --- CDF ---
# Von Mises doesn't have a closed-form CDF. 
# Use Distributions.VonMises for numerical CDF.
function Distributions.cdf(d::VonMisesDistribution, x::Real)
    return Distributions.cdf(Distributions.VonMises(d.μ, d.κ), x)
end

# --- Sampling ---
# Use Distributions.VonMises for sampling   
function Distributions.rand(rng::AbstractRNG, d::VonMisesDistribution)
    return rand(rng, Distributions.VonMises(d.μ, d.κ))
end

# convenience rand when rng::AbstractRNG random generator type is not passed as an argument
Distributions.rand(d::VonMisesDistribution) = rand(Random.GLOBAL_RNG, d)


Statistics.mean(d::VonMisesDistribution) = d.μ
Statistics.var(d::VonMisesDistribution) = 1.0 - (besselj(1, d.κ) / besselj(0, d.κ))  # Calculate the first trigonometric moment R1