module AiryWaves

using Random
using ..SpectralSpreading
using ..AngularSpreading
using ..PhysicalConstants

export AiryState, generate_airy_function

"""
    AiryState
Holds the discrete frequency-direction components of the sea state.
This is the "data" layer.
"""
struct AiryState
    spectrum::DiscreteSpectrum  # Discrete spectrum metadata
    spread::SpreadingModel      # Angular spreading metadata
    nω::Int64                   # Number of frequency bins
    nθ::Int64                   # Number of angle bins
    ω::Vector{Float64}          # Radian frequencies [rad/s]
    k::Vector{Float64}          # Wavenumbers [rad/m]
    θ::Vector{Float64}          # Spreading angles [rad]
    h::Float64                  # Water depth [m]
    seed::Int64                 # Random seed for reproducibility
end

function AiryState(spec::DiscreteSpectrum, spread::AbstractAngularSpreading, h::Real)
    # 1. Radian frequencies from the discrete frequency model
    nω = spec.nbands
    ω_vec = 2π .* get_central_frequencies(spec)
    
    # 2. Angles from the spreading model
    nθ = spread.nθ - 1
    θ_vec = get_central_angles(spread)
    
    # 3. Solve dispersion once for the frequency vector
    k_vec = [solve_wavenumber(w, h) for w in ω_vec]

    return AiryState(spec, spread, nω, nθ, ω_vec, k_vec, θ_vec, Float64(h), Random.seed!())
end


function AiryState(spectrum_model::Symbol, Hs::T, Tp::T,
                   sampling_model::AbstractSampling, fmin, fmax, nf,
                   angular_spreading::Symbol, μ::T, σ::T, θmin::T, θmax::T, nθ::Int,
                   h::T) where {T<:Real}

    # 1. Create continuous spectrum model
    continuous_spectrum = if spectrum_model == :JONSWAP
        JONSWAP(Hs, Tp)
    else
        error("Unsupported spectrum model")
    end

    # 2. Sample spectrum to create discrete spectrum
    spec = DiscreteSpectrum(continuous_spectrum, sampling_model, fmin, fmax, nf)

    # 3. Create the angular spreading model
    spread = SpreadingModel(angular_spreading, μ, σ, θmin, θmax, nθ)

    # 4. Create the AiryState
    return AiryState(spec, spread, h)
end

# --- Seed Management ---

# Helper to handle the "Any" RNG type and seed initialization
function get_seeded_rng(seed::Int64)
    return Random.MersenneTwister(seed) # Or Xoshiro(seed)
end

function change_seed!(state::AiryState, new_seed::Int)
    return AiryState(state.spec, state.spread, state.nω, state.nθ, state.ω, state.k, state.θ, state.h, new_seed)
end

function change_seed!(state::AiryState)
    return change_seed!(state, rand(1:10^9))
end

# -------------------------

# --- INTERNAL DISPERSION SOLVER ---

function solve_wavenumber(ω::Real, h::Real)
    k = ω^2 / g # Deep water guess
    for _ in 1:15
        f = g * k * tanh(k * h) - ω^2
        df = g * tanh(k * h) + g * k * h * (sech(k * h))^2
        dk = f / df
        k -= dk
        if abs(dk) < 1e-8 break end
    end
    return k
end

# --- EVALUATION ---
"""
For each bin (i, j) with center frequency fᵢ and direction θⱼ, we have:
    - Angular Frequency: ωᵢ = 2πfᵢ
    - Wavenumber: kᵢ, solved from the Dispersion Relation: ωᵢ² = g kᵢ tanh(kᵢ h) (where h is the water depth and g is gravity).
    - Amplitude: Aᵢⱼ = √(2 · S(fᵢ) · D(fᵢ, θⱼ) · Δfᵢ · Δθ)
    - Phase: ψᵢⱼ(x, y, z, t) = kᵢ(x cos θⱼ + y sin θⱼ) - ωᵢ t + φᵢⱼ

These elements are combined to compute the sea surface elevation and velocity components at any point (x, y, z) and time t: 
    - Surface Elevation (η) The total displacement of the free surface from the mean water level (z=0):
            η(x, y, t) = Σᵢ Σⱼ Aᵢⱼ cos(ψᵢⱼ)
    - Velocity Potential (Φ) The scalar potential:
            Φ(x, y, z, t) = Σᵢ Σⱼ Aᵢⱼ · g/ωᵢ · [ cosh(kᵢ(z + h)) / cosh(kᵢ h) ] · sin(ψᵢⱼ)

    whose gradient returns the velocity Field (u, v, w):
        · Horizontal (x): 
                u = Σᵢ Σⱼ Aᵢⱼ · ωᵢ · cos(θⱼ) · [ cosh(kᵢ(z + h)) / sinh(kᵢ h) ] · cos(ψᵢⱼ)
        · Horizontal (y): 
                v = Σᵢ Σⱼ Aᵢⱼ · ωᵢ · sin(θⱼ) · [ cosh(kᵢ(z + h)) / sinh(kᵢ h) ] · cos(ψᵢⱼ)
        · Vertical (z): 
                w = Σᵢ Σⱼ Aᵢⱼ · ωᵢ · [ sinh(kᵢ(z + h)) / sinh(kᵢ h) ] · sin(ψᵢⱼ)
"""

"""
    sea_profiles(state::AiryState, x, y, z, t)

The main evaluation engine. Computes η, u, v, w for a single point.
Everything is computed on-the-fly to minimize memory footprint.
"""
function sea_profiles(state::AiryState, x, y, z, t)

    # Amplitudes Matrix  A_ij: (nω × nθ)
    A_ij = get_amplitudes(state)

    # Random phases for all components   ϕ_ij: (nω × nθ)
    ϕ_ij = get_random_phases(state)       

    # Pre-compute phases   ψ_ij: (nω × nθ)
    ψ_ij = state.k .* ( x .* cos.(state.θ) .+ y .* sin.(state.θ) )' .- state.ω .* t .+ ϕ_ij

    # Vertical Coefficients (nω,)
    # We apply the check per-frequency to handle different wave lengths correctly
    kh = state.k .* state.h
    kzh = state.k .* (z + state.h)
    
    # Vectorized conditional: use ifelse. to prevent branching in the loop  (nω,)
    coeff_H = ifelse.(kh .< 20.0, cosh.(kzh) ./ sinh.(kh), exp.(state.k .* z))
    coeff_V = ifelse.(kh .< 20.0, sinh.(kzh) ./ sinh.(kh), exp.(state.k .* z))

    # Compute profiles by summing over all components
    η = sum( A_ij .* cos.(ψ_ij) )
    u = sum( A_ij .* (( state.ω .* coeff_H ) .* cos.(state.θ)') .* cos.(ψ_ij) )
    v = sum( A_ij .* (( state.ω .* coeff_H ) .* sin.(state.θ)') .* cos.(ψ_ij) )
    w = sum( A_ij .* state.ω .* coeff_V .* sin.(ψ_ij) )
    
    return (η=η, u=u, v=v, w=w)
end


function get_amplitudes(state::Airystate)

    # Metadata
    Aω = get_amplitudes(state.spectrum)   # Spectral densities at central frequencies   Aω: (nω,)
    Δθ = get_bandwidth(state.spread)      # Angle bin widths                            Δθ: (nθ,)
    Dθ = get_weights(state.spread)        # Directional spreading weights               Dθ: (nθ,)

    # Amplitudes Matrix  A_ij: (nω × nθ)
    return Aω * sqrt.(Dθ .* Δθ)'

end

function get_random_phases(state::AiryState)
    return 2π .* rand(get_seeded_rng(state.seed), nω, nθ) 
end

end # module