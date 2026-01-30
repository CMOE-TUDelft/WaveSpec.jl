# Bretschneider.jl

"""
    Bretschneider(Hs, Tp)

Creates a Bretschneider spectrum instance (also known as a Modified Pierson-Moskowitz).
    - Hs: Significant wave height (m)
    - Tp: Peak period (s)

The Bretschneider spectrum is self-normalizing and follows an f⁻⁵ slope.
"""
struct Bretschneider{T<:Real} <: AbstractSpectrum
    Hs::T
    Tp::T
    fp::T
end

# --- Primary Constructor ---
function Bretschneider(Hs::Real, Tp::Real)
    T = promote_type(typeof(Hs), typeof(Tp))
    fp = 1.0 / Tp
    return Bretschneider{T}(Hs, Tp, fp)
end

# --- Accessor Functions ---
get_Hs(s::Bretschneider) = s.Hs
get_Tp(s::Bretschneider) = s.Tp


# --- Density Function ---
"""
    get_density(s::Bretschneider, f::Real)

Returns the spectral density S(f) [m²s].
Formula: S(f) = (1.25/4) * (Hs²/fp) * (f/fp)⁻⁵ * exp(-1.25 * (f/fp)⁻⁴)
"""
function get_density(s::Bretschneider, f::Real)
    (f <= 1e-6) && return 0.0
    
    # Relative frequency
    fr = f / s.fp
    
    return 0.3125 * (s.Hs^2 / s.fp)  * fr^-5 * exp(-1.25 * fr^-4)
end