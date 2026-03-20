module AiryWaves

using Random
using ..SpectralSampling
using ..SpectralSpreading
using ..AngularSpreading
using ..PhysicalConstants
using ..ContinuousSpectrums: JONSWAP, RegularWave

export AiryState, generate_sea, get_amplitude, get_random_phases

"""
    AiryState
Holds the discrete frequency-direction components of the sea state.
This is the "data" layer.
"""
struct AiryState
    spectrum::DiscreteSpectralSpreading  # Discrete spectrum metadata
    spread::DiscreteAngularSpreading      # Angular spreading metadata
    nω::Int64                   # Number of frequency bins
    nθ::Int64                   # Number of angle bins
    ω::Vector{Float64}          # Radian frequencies [rad/s]
    k::Vector{Float64}          # Wavenumbers [rad/m]
    θ::Vector{Float64}          # Spreading angles [rad]
    h::Float64                  # Water depth [m]
    seed::Int64                 # Random seed for reproducibility
end

function AiryState(spec::DiscreteSpectralSpreading, spread::DiscreteAngularSpreading, h::Real)
    # 1. Radian frequencies from the discrete frequency model
    nω = spec.nbands
    ω_vec = 2π .* SpectralSpreading.get_central_frequencies(spec)
    
    # 2. Angles from the spreading model
    nθ = spread.nθ - 1
    θ_vec = AngularSpreading.get_central_angles(spread)
    
    # 3. Solve dispersion once for the frequency vector
    k_vec = [solve_wavenumber(w, h) for w in ω_vec]

    return AiryState(spec, spread, nω, nθ, ω_vec, k_vec, θ_vec, Float64(h), abs(rand(Int64)))
end


function AiryState(spec::DiscreteSpectralSpreading, h::Real; θ::Real = 0.0)
    spread = DiscreteAngularSpreading(θ)
    return AiryState(spec, spread, h)
end

function AiryState(spectrum_model::Symbol, Hs::T, Tp::T,
                   sampling_model::AbstractSampling, fmin, fmax, nf,
                   angular_spreading::Symbol, μ::T, σ::T, θmin::T, θmax::T, nθ::Int,
                   h::T) where {T<:Real}

    # 1. Create continuous spectrum model
    continuous_spectrum = if spectrum_model == :JONSWAP
        JONSWAP(Hs, Tp)
    elseif spectrum_model == :RegularWave
        RegularWave(Hs, Tp)
    else
        error("Unsupported spectrum model")
    end

    # 2. Sample spectrum to create discrete spectrum
    spec = DiscreteSpectralSpreading(continuous_spectrum, sampling_model, fmin, fmax, nf)

    # 3. Create the angular spreading model
    spread = DiscreteAngularSpreading(angular_spreading, μ, σ, θmin, θmax, nθ)

    # 4. Create the AiryState
    return AiryState(spec, spread, h)
end

# --- Seed Management ---

# Helper to handle the "Any" RNG type and seed initialization
function get_seeded_rng(seed::Int64)
    return Random.MersenneTwister(seed) # Or Xoshiro(seed)
end

function change_seed!(state::AiryState, new_seed::Int)
    if new_seed < 0 throw(ArgumentError("new_seed must be non-negative (received $new_seed).")) end
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
    generate_sea(state::AiryState, x, y, z, t)

The main evaluation engine. Computes η, u, v, w for a 4D grids (x, y, z, t).
Everything is computed on-the-fly to minimize memory footprint.
The trick to compute the series is to expand all items to 6 dimensional tensors 
with dimensional structure 
                T = (x, y, z, t, ω, θ)
and then contract the last 2 dimensions (ω and θ) to sum the series components.
"""
function generate_sea(state::AiryState, x::AbstractArray{<:R}, y::AbstractArray{<:R}, z::AbstractArray{<:R}, t::AbstractArray{<:R}; vars = [:η, :u, :v, :w]) where {R<:Real}

    # Normalize input to a Vector of Symbols
    requested = vars isa Symbol ? [vars] : vars

    # Error check
    if isempty(requested)
        throw(ArgumentError("No valid variables requested. Choose from :η, :u, :v, :w"))
    end
    
    # Initialize an empty dictionary to store results
    res = Dict{Symbol, Any}()

    # Amplitudes Matrix  A_ij: (nω × nθ)
    A_ij = get_amplitudes(state)

    # Random phases for all components   ϕ_ij: (nω × nθ)
    ϕ_ij = get_random_phases(state)       

    ## Reshape elements to 6 dimensional tensors
    # Reshape Spectral Components to [1, 1, 1, 1, nω, nθ]
    A = reshape(A_ij, 1, 1, 1, 1, state.nω, state.nθ)
    ϕ = reshape(ϕ_ij, 1, 1, 1, 1, state.nω, state.nθ)
    ω = reshape(state.ω, 1, 1, 1, 1, state.nω, 1)
    θ = reshape(state.θ, 1, 1, 1, 1, 1, state.nθ)
    k = reshape(state.k, 1, 1, 1, 1, state.nω, 1)

    # Reshape Evaluation Coordinates
    X = reshape(x, :, 1, 1, 1, 1, 1) # Dim 1
    Y = reshape(y, 1, :, 1, 1, 1, 1) # Dim 2
    Z = reshape(z, 1, 1, :, 1, 1, 1) # Dim 3
    T = reshape(t, 1, 1, 1, :, 1, 1) # Dim 4

    # Phase Tensor Construction (nx × ny × nz × nt × nω × nθ)
    # ψ = k*(x*cosθ + y*sinθ) - ωt + ϕ
    ψ = k .* ( X .* cos.(θ) .+ Y .* sin.(θ) ) .- (ω .* T) .+ ϕ

    # Calculations (Independent 'if' blocks so all requested vars are computed) -> Compute profiles by summing over all components
    if :η in requested
        res[:η] = dropdims(sum( A .* cos.(ψ), dims=(5,6)), dims=(5,6))
    end
    
    if :u in requested || :v in requested || :w in requested
        # Vertical Coefficients (Expanded to handle Z and k simultaneously)
        kh = state.k .* state.h  # (nω, )
        # Reshape k and kh for Dim 5 (ω)
        k_reshaped  = reshape(state.k, 1, 1, 1, 1, state.nω, 1)
        kh_reshaped = reshape(kh, 1, 1, 1, 1, state.nω, 1)
        h = state.h 
    end

    if :u in requested || :v in requested 
        # Precompute coefficients only if  needed: coeff_H will be (1 × 1 × nz × 1 × nω × 1)
        coeff_H = ifelse.(kh_reshaped .< 20.0, 
                    cosh.(k_reshaped .* (Z .+ h)) ./ sinh.(kh_reshaped), 
                    exp.(k_reshaped .* Z))

        if :u in requested
            # Compute profile by summing over all components
            res[:u] = dropdims(sum( A .* (ω .* coeff_H .* cos.(θ)) .* cos.(ψ), dims=(5, 6)), dims=(5, 6))
        end

        if :v in requested
            # Compute profile by summing over all components
            res[:v] = dropdims(sum( A .* (ω .* coeff_H .* sin.(θ)) .* cos.(ψ), dims=(5, 6)), dims=(5, 6))
        end
    end

    if :w in requested
        # Precompute coefficients only if needed: coeff_V will be (1 × 1 × nz × 1 × nω × 1)
        # We apply the check per-frequency to handle different wave lengths correctly
        coeff_V = ifelse.(kh_reshaped .< 20.0, 
            sinh.(k_reshaped .* (Z .+ h)) ./ sinh.(kh_reshaped), 
            exp.(k_reshaped .* Z))
        # Compute profile by summing over all components
        res[:w] = dropdims(sum( A .* (ω .* coeff_V) .* sin.(ψ), dims=(5, 6)), dims=(5, 6))
    end

    return (; res...)
end


function get_amplitudes(state::AiryState)

    # Metadata
    Aω = SpectralSpreading.get_amplitudes(state.spectrum)   # Spectral densities at central frequencies   Aω: (nω,)
    Δθ = AngularSpreading.get_bandwidths(state.spread)     # Angle bin widths                            Δθ: (nθ,)
    Dθ = AngularSpreading.get_weights(state.spread)        # Directional spreading weights               Dθ: (nθ,)

    # Amplitudes Matrix  A_ij: (nω × nθ)
    return Aω * sqrt.(Dθ .* Δθ)'

end

function get_random_phases(state::AiryState)
    return 2π .* rand(get_seeded_rng(state.seed), state.nω, state.nθ) 
end

end # module