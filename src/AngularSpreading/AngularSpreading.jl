module AngularSpreading

using Distributions
using Random
using ..Truncation

include("CosinePower.jl")
include("VonMises.jl")
include("DonelanBanner.jl")

export SpreadingModel, get_seeded_rng
export get_angle, get_angles, 
       get_central_angle, get_central_angles, 
       get_bandwidth, get_bandwidths, get_weights

"""
    SpreadingModel{D<:UnivariateDistribution, T<:Real}

A structure to define angular spreading models for wave spectra.
- `distribution`: A truncated univariate distribution defining the angular spreading.
- `nθ`: Number of discrete angles to sample from the distribution.
"""
struct SpreadingModel{D<:UnivariateDistribution, T<:Real}
    distribution::TruncatedModel{D, T}  # Truncated angular distribution
    nθ::Int64                           # Number of discrete angles
    seed::Int64                         # Random seed for reproducibility
end

"""
    SpreadingModel(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int) 

Creates a SpreadingModel by truncating the given univariate distribution between [a, b].
"""
function SpreadingModel(dist::UnivariateDistribution, a::Real, b::Real, nθ::Int)
    truncated = TruncatedModel(dist, a, b)
    return SpreadingModel{typeof(dist), Int64, Int64}(truncated, nθ, Random.seed!())
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
    
    return SpreadingModel{typeof(base_dist), Int64, Int64}(truncated_dist, nθ, Random.seed!())
end

# --- Seed Management ---

# Helper to handle the "Any" RNG type and seed initialization
function get_seeded_rng(seed::Int64)
    return Random.MersenneTwister(seed) # Or Xoshiro(seed)
end

function change_seed!(sm::SpreadingModel, new_seed::Int)
    return SpreadingModel(sm.distribution, sm.nθ, new_seed)
end

function change_seed!(sm::SpreadingModel)
    return change_seed!(sm, rand(1:10^9))
end

# -----------------------

# --- Attributes & Sampling ---

"""
    get_angles(sm::SpreadingModel, rng=AbstractRNG)
Returns n_theta random angles sampled from the truncated distribution.
"""
function get_angles(sm::SpreadingModel)
    return sort(rand(get_seeded_rng(sm.seed), sm.distribution, sm.nθ))
end

function get_angles(sm::SpreadingModel, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r), sm.nθ)
    (start_idx > end_idx) && return Float64[]
    return get_angles(sm)[start_idx:end_idx]
end

function get_angle(sm::SpreadingModel, idx::Int)
    (idx < 1 || idx > sm.nθ) && Float64[]
    return get_angles(sm)[idx]
end

"""
    get_central_angles(sm::SpreadingModel)
Returns the nθ - 1 central angles between the sampled edges.
"""
function get_central_angles(sm::SpreadingModel)
    angles = get_angles(sm)
    # Compute midpoints: (θ_{j} + θ_{j+1}) / 2
    return (angles[1:end-1] .+ angles[2:end]) ./ 2.0
end

function get_central_angles(sm::SpreadingModel, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_central_angles(sm)[start_idx:end_idx]
end

function get_central_angle(sm::SpreadingModel, idx::Int)
    (idx < 1 || idx > sm.nθ - 1) && return Float64[]
    return get_central_angles(sm)[idx]
end

"""
    get_bandwidths(sm::SpreadingModel)
Returns the Δθ for the nθ - 1 bins.
"""
function get_bandwidths(sm::SpreadingModel)
    return diff(get_angles(sm))
end

function get_bandwidths(sm::SpreadingModel, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_bandwidths(sm)[start_idx:end_idx]
end

function get_bandwidth(sm::SpreadingModel, idx::Int)
    (idx < 1 || idx > sm.nθ - 1) && return Float64[]
    return get_bandwidths(sm)[idx]
end

"""
    get_weights(sm::SpreadingModel)

Returns the PDF values at central bins angles, corrected for discrete sampling.
"""
function get_weights(sm::SpreadingModel)
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

function get_weights(sm::SpreadingModel, r::AbstractRange)
    start_idx = max(1, first(r))
    end_idx   = min(last(r)-1, sm.nθ-1)
    (start_idx > end_idx) && return Float64[]
    return get_weights(sm)[start_idx:end_idx]
end

end # module AngularSpreading