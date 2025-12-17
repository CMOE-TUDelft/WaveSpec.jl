module AngularSpreading

using Distributions
using Random
using ..Truncation

include("CosinePower.jl")
include("VonMises.jl")
include("DonelanBanner.jl")

export SpreadingModel, get_angles, get_weights, get_angles_weights

"""
    SpreadingModel{D<:UnivariateDistribution, T<:Real}

A structure to define angular spreading models for wave spectra.
- `distribution`: A truncated univariate distribution defining the angular spreading.
- `nθ`: Number of discrete angles to sample from the distribution.
"""
struct SpreadingModel{D<:UnivariateDistribution, T<:Real}
    distribution::TruncatedModel{D, T}
    nθ::Int
end

"""
    SpreadingModel(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int) 

Creates a SpreadingModel by truncating the given univariate distribution between [a, b].
"""
function SpreadingModel(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int)
    truncated = TruncatedModel(dist, a, b)
    return SpreadingModel{typeof(dist), typeof(a)}(truncated, nθ)
end

"""
    SpreadingModel(model_type::Symbol, μ, σ, a, b, nθ)

Factory constructor.
- :uniform -> Uniform 
- :normal  -> Gaussian 
- :cosinepow  -> Cosine Power (2s).
- :vonmises -> Circular Normal.
- :donelan -> Donelan-Banner (sech²).
"""
function SpreadingModel(model_type::Symbol, μ::T, σ::T, a::T, b::T, nθ::Int) where {T<:Real}
    # μ: mean direction (radians)
    # σ: standard deviation (radians)
    # a, b: truncation limits (radians)
    # nθ: number of discrete angles to sample
    
    base_dist = if model_type == :uniform
        Uniform(a, b)
    
    elseif model_type == :normal
        Normal(μ, σ)
    
    elseif model_type == :cosinepow
        # For Cosine-2s, the variance relates to n (or 2s). 
        # A common approximation: n ≈ 2/(σ^2) - 2
        n_val = max(0.1, 2.0 / (σ^2) - 2.0)
        CosinePowerDistribution(μ, n_val)
    
    elseif model_type == :vonmises
        # Unification: σ² ≈ 1/κ  => κ = 1/σ²
        κ_val = 1.0 / (σ^2)
        VonMisesDistribution(μ, κ_val)
    
    elseif model_type == :donelan
        # Unification: σ² ≈ π²/(3β²) => β = π / (sqrt(3) * σ)
        β_val = π / (sqrt(3.0) * σ)
        DonelanBannerDistribution(μ, β_val)
    
    else
        error("Model type :$model_type not recognized.")
    end

    if model_type == :uniform
        # Uniform distribution does not need truncation
        truncated_dist = base_dist
    else
        # Apply truncation for other distributions
        truncated_dist = TruncatedModel(base_dist, a, b)
    end
    
    return SpreadingModel{typeof(base_dist), T}(truncated_dist, nθ)
end

# --- Attributes & Sampling ---

"""
    get_angles(sm::SpreadingModel; rng=Random.GLOBAL_RNG)
Returns n_theta random angles sampled from the truncated distribution.
"""
function get_angles(sm::SpreadingModel; rng::AbstractRNG=Random.GLOBAL_RNG)
    # Sample nθ angles from the truncated distribution
    θⱼ = sort([rand(rng, sm.distribution) for _ in 1:sm.nθ])
    # Compute bins' widths
    Δθⱼ = diff(θⱼ) 
    push!(Δθⱼ,Δθⱼ[end]) # last bin = second last bin

    return θⱼ, Δθⱼ
end

"""
    get_weights(sm::SpreadingModel, angles::Vector{Real})
Returns the PDF values at specific angles. Because the angles are sampled 
from the distribution itself, the weights for a Monte Carlo sum would 
actually be 1/n_theta, but we provide the PDF values here for completeness.
"""
function get_weights(sm::SpreadingModel, θⱼ::Vector{<:Real})
    # PDF values at sampled angles
    weights = [pdf(sm.distribution, θ) for θ in θⱼ]
    # Compute bins' widths
    Δθⱼ = diff(θⱼ) 
    push!(Δθⱼ,Δθⱼ[end]) # last bin = second last bin
    # Compute total discrete probability (check normalisation) 
    Σ_θⱼ = sum(weights.*Δθⱼ)
    
    return Δθⱼ, weights/Σ_θⱼ
end


"""
    get_angles_weights(sm::SpreadingModel; rng=Random.GLOBAL_RNG)

Returns a Tuple of (angles, weights).
- `angles`: n_theta samples drawn from the distribution (Importance Sampling).
- `weights`: The PDF values at those specific samples.

This is the primary function for 3D Spectral initialization.
"""
function get_angles_weights(sm::SpreadingModel; rng::AbstractRNG=Random.GLOBAL_RNG)
    # 1. Sample angles based on the distribution shape
    θⱼ, Δθⱼ = get_angles(sm; rng=rng)
    # 2. Calculate the PDF weights at these specific points
    weights = [pdf(sm.distribution, θ) for θ in θⱼ]
    # Compute total discrete probability (check normalisation) 
    Σ_θⱼ = sum(weights.*Δθⱼ)

    return θⱼ, Δθⱼ, weights/Σ_θⱼ
end


end # module AngularSpreading