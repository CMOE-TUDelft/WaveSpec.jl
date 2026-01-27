# Donelan.jl

"""
    Donelan(Hs, Tp, U10)

Creates a Donelan spectrum instance.
    - Hs: Significant wave height (m)
    - Tp: Peak period (s)
    - U10: Wind speed at 10m height (m/s)

The Donelan spectrum follows an f⁻⁴ slope in the high-frequency equilibrium range.
"""
struct Donelan{T<:Real} <: AbstractSpectrum
    Hs::T
    Tp::T
    U10::T
    fp::T
    γ::T    # Peak enhancement factor
    Ag::T   # Normalization factor
    σ₁::T
    σ₂::T
end

# --- Primary Constructor ---
function Donelan(Hs::Real, Tp::Real, U10::Real)
    T = promote_type(typeof(Hs), typeof(Tp), typeof(U10))
    fp = 1.0 / Tp
    
    # Donelan's specific γ estimation based on wave age
    # cp is the phase speed at the peak frequency
    cp = 9.81 / (2π * fp)
    U_cp = U10 / cp  # Inverse wave age
    
    # Peak enhancement γ calculation (Donelan et al., 1985)
    if U_cp < 1.0
        γ = 1.7
    else
        γ = 1.7 + 6.0 * log10(U_cp)
    end
    γ = clamp(γ, 1.7, 7.0)

    # Standard Donelan width parameters
    σ₁ = 0.08
    σ₂ = 0.10

    # Like JONSWAP, this requires numerical normalization to ensure m0 = (Hs/4)²
    Ag = compute_donelan_normalization(γ, fp)

    return Donelan{T}(Hs, Tp, U10, fp, γ, Ag, σ₁, σ₂)
end

# --- 2. Accessor Functions ---
get_Hs(s::Donelan) = s.Hs
get_Tp(s::Donelan) = s.Tp

"""
    get_density(s::Donelan, f::Real)

Returns the spectral density S(f). Note the f⁻⁴ dependency.
"""
function get_density(s::Donelan, f::Real)
    (f <= 1e-6) && return 0.0
    
    fr = f / s.fp
    σ = (fr <= 1.0) ? s.σ₁ : s.σ₂
    b = exp(-0.5 * ((fr - 1.0) / σ)^2)
    
    # Donelan Formulation: S(f) ∝ f⁻⁴ * exp(-(fp/f)²) * γ^b
    # The exponent for the decay is -2.0, and the slope is fr^-4
    return s.Ag * ((s.Hs / 4.0)^2 / s.fp) * (fr^-4 * exp(-fr^-2)) * s.γ^b
end

"""
    compute_donelan_normalization(γ, fp)

Numerical normalization for the Donelan shape.
"""
function compute_donelan_normalization(γ::Real, fp::Real)
    function donelan_shape(f::Real)
        (f <= 1e-6) && return 0.0
        fr = f / fp
        σ = (fr <= 1.0) ? 0.08 : 0.10
        b = exp(-0.5 * ((fr - 1.0) / σ)^2)
        # Unit shape (Hs=4, fp=1)
        return (1.0/fp) * (fr^-4 * exp(-fr^-2)) * γ^b
    end

    # Integrate from near 0 to high frequency
    total_area = IntegrateGaussQuad(donelan_shape, order=4, a=1e-4, b=15.0, n=250)
    return 1.0 / total_area
end