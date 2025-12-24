# JONSWAP.jl

"""
    JONSWAP(Hs, Tp, γ=nothing)

Creates a JONSWAP (JOint North Sea WAve Project) spectrum instance with given:
    - significant wave height Hs (m) -> Average height of the highest third of waves
    - peak period Tp (s) -> Period of waves with the most energy
    - optional peak enhancement factor γ -> Controls the sharpness/height of the spectrum's peak,
If γ is not provided, estimated based on Hs and Tp.
"""


struct JONSWAP{T<:Real} <: AbstractSpectrum
    Hs::T   # Significant wave height
    Tp::T   # Peak period
    fp::T   # Peak frequency
    γ::T    # Peakedness factor     # γ = 3.3   -> 5.870*fpeak^0.86  (See holthuis pg.162)     
    Ag::T   # Normalization factor
    ωp::T   # Peak angular frequency
    σ₁::T   # Spectral width parameter for frequencies below peak
    σ₂::T   # Spectral width parameter for frequencies above peak
end

# --- 1. Full Constructor (All parameters provided) ---
# This is your primary entry point
function JONSWAP(Hs::T, Tp::T, γ::T) where {T<:Real}
    fp = 1.0 / Tp                       
    Ag = compute_normalization(γ)
    ωp = 2π / Tp
    σ₁ = 0.07           # 0.0547*fpeak^0.32
    σ₂ = 0.09           # 0.0783*fpeak^0.16
    return JONSWAP{T}(Hs, Tp, fp, γ, Ag, ωp, σ₁, σ₂)
end

# --- 2. Partial Constructor (Hs and Tp only) ---
# Automatically estimates γ
function JONSWAP(Hs::Real, Tp::Real)
    est_γ = estimate_γ(Hs, Tp)
    return JONSWAP(Hs, Tp, est_γ)
end

# --- 3. Partial Constructor (Tp and γ only) ---
# Automatically estimates Hs using your estimate_Hs utility
# Useful for "Standard" sea states where height is defined by period
function JONSWAP(;Tp::Real, γ::Real=3.3)
    est_Hs = estimate_Hs(Tp, γ)
    return JONSWAP(est_Hs, Tp, γ)
end

# --- 4. Partial Constructor (Hs and γonly) ---
# Assuming a standard steepness to estimate a reasonable Tp, then estimating γ   
function JONSWAP(Hs::Real, γ::Real=3.3)
    est_Tp = estimate_Tp(Hs, γ)
    return JONSWAP(Hs, est_Tp, γ)
end

"""
    estimate_γ(Hs, Tp)

Estimation of γ, given Hs and Tp
"""
function estimate_γ(Hs::Real, Tp::Real)
    x = Tp / sqrt(Hs)
    if x < 5.143
        D = 0.036 - 0.0056 * x
        γ = exp(3.484 * (1.0 - 0.1975 * D * x^4))
    else
        γ = 1.0
    end
    return clamp(γ, 1.0, 7.0)
end

"""
    estimate_Hs(Tp, γ=3.3)

Estimation of Hs, given Tp and gamma
"""
function estimate_Hs(Tp::Real, γ = 3.3)
  return (0.11661 + 0.01581*γ - 0.00065*γ*γ)*Tp*Tp
end

"""
    estimate_Tp(Hs, γ=3.3)
Estimation of Tp, given Hs and gamma
"""
function estimate_Tp(Hs::Real, γ=3.3)
    a = 0.11661 + 0.01581*γ - 0.00065*γ*γ
    return sqrt(Hs / a)
end


"""
    get_density(s::JONSWAP, f::Real)

Returns the spectral density S(f) in [m²s].
Assumes s.Ag is the normalization factor for the f-based JONSWAP.
"""
function get_density(s::JONSWAP, f::Real)
    # 1. Guard against zero or negative frequency
    (f <= 1e-6) && return 0.0
    # Compute relative frequency (dimensionless)
    fr = f / s.fp        
    # Shape parameters
    b = get_γ_exponent(s, fr)
    # Spectral density (energy-normalized form)
    return s.Ag * ((s.Hs / 4)^2 / s.fp) * (5.0 * fr^-5 * exp(-1.25 * fr^-4)) * s.γ^b
end


"""
    get_γ_exponent(s::JONSWAP, fr::Real)

Returns the γ exponent b for a given relative frequency fr = f/fp.
"""
function get_γ_exponent(s::JONSWAP, fr::Real) 
    σ = fr <= 1.0 ? s.σ₁ : s.σ₂
    b = exp(-0.5 * ((fr - 1.0) / σ)^2)
    return b
end

"""
    get_alpha(s::JONSWAP)

Estimation of α, given JONSWAP parameters
"""
function get_alpha(s::JONSWAP)
    return s.Hs^2 * s.fp^4 * s.Ag * π^4 / g^2
end


"""
    compute_jonswap_norm(γ)

Calculates the normalization factor Ag needed to ensure m₀ = (Hs/4)².
The most accurate way is to integrate the dimensionless shape function. 
Because fr = f/fp, the peak always occurs at fr = 1.0. This makes the integration independent of the actual H_s or T_p values.
This way, we don't have to worry about the specific frequency range of a given sea state (e.g., Tp=5s vs Tp=15s) because the "shape" is always centered at 1.0.
"""
# Compute normalization factor Ag
function compute_normalization(γ::Real)
    # 1. Define the dimensionless "unit" JONSWAP spectrum (Hs=4 and fp=1)
    # G(fr) = 5 * fr^-5 * exp[ -(5/4) * fr^-4 ] * γ^{ exp[ -0.5 * ((fr - 1)/σ)^2 ] }
    # capturing 'γ' from the outer scope
    function G_fr(fr::Real)
        (fr <= 1e-6) && return 0.0
        # Standard JONSWAP constants
        σ = (fr <= 1.0) ? 0.07 : 0.09
        # Peak enhancement exponent
        b = exp(-0.5 * ((fr - 1.0) / σ)^2)
        # Dimensionless JONSWAP shape
        return (5.0 * fr^-5 * exp(-1.25 * fr^-4)) * γ^b
    end

    # Integrate the "unit" JONSWAP spectrum from ~0 to a very high fr to capture the tail
    # Lower integration bound   a = 1e-4 * fp = 1e-4 * 1.0 = 1e-4
    # Upper integration bound   b = 20.0 * fp = 20.0 * 1.0 = 20.0
    total_area = gaussQuad1D_4(f=G_fr, a=1e-4, b=20.0, N=200)

    return 1.0 / total_area
end

