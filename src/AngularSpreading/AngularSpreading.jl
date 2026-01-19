module AngularSpreading

using Distributions
using Random
using RecipesBase
using Statistics
using SpecialFunctions
using ..Truncation

include("CosinePower.jl")
include("VonMises.jl")
include("DonelanBanner.jl")

export DiscreteAngularSpreading, get_seeded_rng
export get_angle, get_angles, 
       get_central_angle, get_central_angles, 
       get_bandwidth, get_bandwidths, get_weights
export plot_model

"""
    DiscreteAngularSpreading{D<:UnivariateDistribution, T<:Real}

A structure to define angular spreading models for wave spectra.
- `distribution`: A truncated univariate distribution defining the angular spreading.
- `nθ`: Number of discrete angles to sample from the distribution.
"""
struct DiscreteAngularSpreading{D<:UnivariateDistribution, T<:Real}
    distribution::TruncatedModel{D, T}  # Truncated angular distribution
    nθ::Int64                           # Number of discrete angles
    seed::Int64                         # Random seed for reproducibility
end

"""
    DiscreteAngularSpreading(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int) 

Creates a DiscreteAngularSpreading by truncating the given univariate distribution between [a, b].
"""
function DiscreteAngularSpreading(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int; units::Symbol=:radians)
    if units == :degrees
        a = a * (π/180)
        b = b * (π/180)
    end
    truncated = TruncatedModel(dist, a, b)
    return DiscreteAngularSpreading(truncated, nθ, abs(rand(Int64)))
end

"""
    DiscreteAngularSpreading(model_type::Symbol, μ, σ, a, b, nθ)

Factory constructor.
- :uniform -> Uniform 
- :normal  -> Gaussian 
- :cosinepow  -> Cosine Power (2s).
- :vonmises -> Circular Normal.
- :donelan -> Donelan-Banner (sech²).
"""
function DiscreteAngularSpreading(model_type::Symbol, μ::T, σ::T, a::T, b::T, nθ::Int; units::Symbol=:radians) where {T<:Real}
    # μ: mean direction [rad]
    # σ: standard deviation [rad]
    # a, b: truncation limits [rad]
    # nθ: number of discrete angles to sample

    if units == :degrees
        μ = μ * (π/180)
        σ = σ * (π/180)
        a = a * (π/180)
        b = b * (π/180)
    end
    
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

    # Apply truncation
    truncated_dist = TruncatedModel(base_dist, a, b)
    return DiscreteAngularSpreading(truncated_dist, nθ, abs(rand(Int64)))
end

function DiscreteAngularSpreading(dist::TruncatedModel{D, T}, nθ::Integer, seed::Integer) where {D, T}
    return DiscreteAngularSpreading{D, T}(dist, Int64(nθ), Int64(seed))
end


"""
    DiscreteAngularSpreading(θ)

Default constructor for a long-crested sea (no angular spreading) in θ direction.
Sets nθ = 2 (2 samples to define a single bin) and uses a narrow Uniform distribution centered at θ.
"""
function DiscreteAngularSpreading(θ::Real; seed::Int64 = abs(rand(Int64)), units::Symbol=:radians)
    if units == :degrees
        θ = θ * (π/180)
    end
    # Define a dummy narrow distribution with a single sample θ 
    # We use a very small range so that the 'mean' is effectively 0
    ϵ = 0.001
    dist = Uniform(θ  - ϵ, θ + ϵ)
    truncated_dist = TruncatedModel(dist, -π, π)
    
    # We force nθ = 1 for long-crested waves
    return DiscreteAngularSpreading(truncated_dist, 2, seed)
end

function DiscreteAngularSpreading()
    return DiscreteAngularSpreading(0.0)
end


# --- Seed Management ---

# Helper to handle the "Any" RNG type and seed initialization
function get_seeded_rng(seed::Int64)
    return Random.MersenneTwister(seed) # Or Xoshiro(seed)
end

function change_seed!(sm::DiscreteAngularSpreading, new_seed::Int)
    if new_seed < 0 throw(ArgumentError("new_seed must be non-negative (received $new_seed).")) end
    return DiscreteAngularSpreading(sm.distribution, sm.nθ, new_seed)
end

function change_seed!(sm::DiscreteAngularSpreading)
    return change_seed!(sm, rand(1:10^9))
end

# -----------------------

# --- Attributes & Sampling ---

"""
    get_angles(sm::DiscreteAngularSpreading, rng=AbstractRNG)
Returns n_theta random angles sampled from the truncated distribution.
"""
function get_angles(sm::DiscreteAngularSpreading)
    return sort(rand(get_seeded_rng(sm.seed), sm.distribution, sm.nθ))
end

function get_angles(sm::DiscreteAngularSpreading, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r), sm.nθ)
    (start_idx > end_idx) && return Float64[]
    return get_angles(sm)[start_idx:end_idx]
end

function get_angle(sm::DiscreteAngularSpreading, idx::Int)
    (idx < 1 || idx > sm.nθ) && Float64[]
    return get_angles(sm)[idx]
end

"""
    get_central_angles(sm::DiscreteAngularSpreading)
Returns the nθ - 1 central angles between the sampled edges.
"""
function get_central_angles(sm::DiscreteAngularSpreading)
    angles = get_angles(sm)
    # Compute midpoints: (θ_{j} + θ_{j+1}) / 2
    return (angles[1:end-1] .+ angles[2:end]) ./ 2.0
end

function get_central_angles(sm::DiscreteAngularSpreading, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_central_angles(sm)[start_idx:end_idx]
end

function get_central_angle(sm::DiscreteAngularSpreading, idx::Int)
    (idx < 1 || idx > sm.nθ - 1) && return Float64[]
    return get_central_angles(sm)[idx]
end

"""
    get_bandwidths(sm::DiscreteAngularSpreading)
Returns the Δθ for the nθ - 1 bins.
"""
function get_bandwidths(sm::DiscreteAngularSpreading)
    return diff(get_angles(sm))
end

function get_bandwidths(sm::DiscreteAngularSpreading, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_bandwidths(sm)[start_idx:end_idx]
end

function get_bandwidth(sm::DiscreteAngularSpreading, idx::Int)
    (idx < 1 || idx > sm.nθ - 1) && return Float64[]
    return get_bandwidths(sm)[idx]
end

"""
    get_weights(sm::DiscreteAngularSpreading)

Returns the PDF values at central bins angles, corrected for discrete sampling.
"""
function get_weights(sm::DiscreteAngularSpreading)
    # Compute central angles
    θⱼ = get_central_angles(sm)
    # PDF values at central angles
    weights = [pdf(sm.distribution, θ) for θ in θⱼ]
    # Compute bins' widths
    Δθⱼ = get_bandwidths(sm)
    # Compute total discrete probability (check normalisation) 
    Σ_θⱼ = sum(weights.*Δθⱼ)

    return weights/Σ_θⱼ
end

function get_weights(sm::DiscreteAngularSpreading, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_weights(sm)[start_idx:end_idx]
end


# --- Plots and Visualisation ---
# Plot using light package RecipesBase:
#   1. import Plots in script
#   2. Call this function as plot(DiscreteAngularSpreading)

@recipe function f(sm::DiscreteAngularSpreading; n_points=200)
    # Plot Attributes
    title  := "Angular spreading discretization"
    xlabel := "Angle θ (rads)"
    ylabel := "Probability density"
    grid   := false

    # Get plot range
    a, b = sm.distribution.a, sm.distribution.b
    θ_range = range(a, b, length=n_points)
    
    # Define the first series: The Continuous Line
    @series begin
        label := "Continuous PDF"
        seriestype := :path
        fillrange := 0
        fillalpha := 0.2
        linecolor := :black
        θ_range, [pdf(sm.distribution, θ) for θ in θ_range]
    end

    # Define the second series: The Stems/Samples
    θ_samples = get_angles(sm)
    weights = [pdf(sm.distribution, θ) for θ in θ_samples]
    
    @series begin
        label := "Discrete Samples"
        seriestype := :scatter
        marker := :circle
        markersize := 3
        markercolor := :red
        θ_samples, weights
    end
    
    # Stem lines logic
    for (t, w) in zip(θ_samples, weights)
        @series begin
            label := false
            seriestype := :path
            linecolor := :red
            [t, t], [0, w]
        end
    end
end

end # module AngularSpreading