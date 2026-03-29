module AiryWaves

using Random
using ..SpectralSampling
using ..SpectralSpreading
using ..AngularSpreading
using ..PhysicalConstants
using Interpolations
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
        error("Unsupported constructor for RegularWave. Please use AiryState(RegularWave(H, T), h) instead.")
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

"""
    generate_interpolable_sea(state::AiryState, x, y, z, t; vars=[:η, :u, :v, :w], interp=:linear)

Computes the sea fields on the provided regular grids and returns a NamedTuple
of interpolation functions for each requested variable. Each interpolant is a
callable taking `(x, y, z, t)` and returning the interpolated value.

Arguments
- `state`: an `AiryState` produced by `AiryState(...)` or constructed manually.
- `x, y, z, t`: regular 1D coordinate vectors used to compute the fields.
- `vars`: which fields to compute, default `[:η, :u, :v, :w]`.
- `interp`: interpolation kernel, `:linear` (default), `:cubic`, or `:nearest`.

Examples

```julia
using WaveSpec

# (1) Construct an AiryState (example uses a regular monochromatic wave)
spec = WaveSpec.ContinuousSpectrums.RegularWave(1.0, 5.0)                      # H=1 m, T=5 s
ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec)                # discrete spectrum
spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)              # unidirectional
as = WaveSpec.AiryWaves.AiryState(ds, spread, 50.0)                          # water depth 50 m

# (2) Choose evaluation grid
x = range(0.0, stop=10.0, length=11)
y = [0.0]
z = [0.0]
t = 0.0:0.1:10.0

# (3) Get interpolants (linear by default)
itps = WaveSpec.AiryWaves.generate_interpolable_sea(as, collect(x), collect(y), collect(z), collect(t); vars=[:η])

# (4) Evaluate η at an arbitrary point (not necessarily on the grid)
η_val = itps[:η](1.23, 0.0, 0.0, 2.57)
```

"""
function generate_interpolable_sea(state::AiryState, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, z::AbstractVector{<:Real}, t::AbstractVector{<:Real}; vars = [:η, :u, :v, :w], interp::Symbol = :linear)
    # Compute sea on the provided grid
    res = generate_sea(state, x, y, z, t; vars = vars)
    interp_dict = Dict{Symbol, Function}()
    axes = (x, y, z, t)
    for v in vars
        if v in keys(res)
            vals = res[v]
            # Ensure the computed array has the same dimensionality as the axes tuple
            target_lens = (length(x), length(y), length(z), length(t))
            cur_shape = size(vals)
            # If the value array has fewer dims, add singleton dimensions at the end
            if length(cur_shape) < 4
                newshape = Tuple(vcat(collect(cur_shape), ones(Int, 4 - length(cur_shape))))
                vals = reshape(vals, newshape)
                cur_shape = size(vals)
            end
            # Repeat singleton dimensions to match the requested axes lengths
            reps = ntuple(i -> cur_shape[i] == target_lens[i] ? 1 : target_lens[i], 4)
            if any(r -> r != 1, reps)
                vals = repeat(vals, reps...)
            end
            # Support several kernels via Interpolations.jl
            if interp == :linear
                itp = Interpolations.interpolate((x, y, z, t), vals, Interpolations.Gridded(Interpolations.Linear()))
            elseif interp == :nearest
                # Build a simple nearest-neighbor wrapper without Interpolations.jl
                itp = nothing
                local_vals = vals
                function nearest_wrapper(xq::Real, yq::Real, zq::Real, tq::Real)
                    # find nearest indices along each axis
                    ix = findmin(abs.(x .- xq))[2]
                    iy = findmin(abs.(y .- yq))[2]
                    iz = findmin(abs.(z .- zq))[2]
                    it = findmin(abs.(t .- tq))[2]
                    return local_vals[ix, iy, iz, it]
                end
                interp_dict[v] = nearest_wrapper
                continue
            else
                itp = nothing
            end

            if itp !== nothing
                function wrapped(xq::Real, yq::Real, zq::Real, tq::Real)
                    xc = clamp(xq, first(x), last(x))
                    yc = clamp(yq, first(y), last(y))
                    zc = clamp(zq, first(z), last(z))
                    tc = clamp(tq, first(t), last(t))
                    return itp(xc, yc, zc, tc)
                end
                interp_dict[v] = wrapped
            else
                error("Unsupported interpolation method: $interp")
            end
        end
    end
    return (; interp_dict...)
end

end # module