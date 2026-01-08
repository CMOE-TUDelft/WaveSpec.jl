# OchiHubble.jl
using SpecialFunctions: gamma

"""
    OchiHubble(Hs1, Tp1, λ1, Hs2, Tp2, λ2)

Bimodal spectrum defined by two independent Ochi-Hubble components.
Typically Component 1 is Swell and Component 2 is Wind Sea.
"""
struct OchiHubble{T<:Real} <: AbstractSpectrum
    Hs1::T          # Lower frequency significant wave height
    Tp1::T          # Lower frequency peak period
    λ1::T           # Lower frequency peak factor
    fp1::T          # Lower frequency peak frequency
    Hs2::T          # Higher frequency significant wave height
    Tp2::T          # Higher frequency peak period
    λ2::T           # Higher frequency peak factor
    fp2::T          # Higher frequency peak frequency
    total_Hs::T     # Total significant wave height
end

# --- 1. Primary Constructor ---
function OchiHubble(Hs1::Real, Tp1::Real, λ1::Real, 
                    Hs2::Real, Tp2::Real, λ2::Real)
                    
    T = promote_type(typeof(Hs1), typeof(Tp1), typeof(λ1), 
                     typeof(Hs2), typeof(Tp2), typeof(λ2))
    
    fp1 = 1.0 / Tp1
    fp2 = 1.0 / Tp2
    total_Hs = sqrt(Hs1^2 + Hs2^2)
    
    return OchiHubble{T}(Hs1, Tp1, λ1, Hs2, Tp2, λ2, fp1, fp2, total_Hs)
end

# --- 2. Keyword Constructor ---
function OchiHubble(;Hs1=1.0, Tp1=14.0, λ1=1.0, 
                     Hs2=2.0, Tp2=7.0, λ2=2.0)
    return OchiHubble(Hs1, Tp1, λ1, Hs2, Tp2, λ2)
end

"""
    get_density(s::OchiHubble, f::Real)

Returns the spectral density S(f) [m²s] by summing the two components.
"""
function get_density(s::OchiHubble, f::Real)
    (f <= 1e-6) && return 0.0
    
    # Component 1
    S1 = compute_ochi_component(f, s.Hs1, s.fp1, s.λ1)
    # Component 2
    S2 = compute_ochi_component(f, s.Hs2, s.fp2, s.λ2)
    
    return S1 + S2
end

"""
    compute_ochi_component(f, Hs, fp, λ)

Calculates the f-based density for a single Ochi component.
Formula: S(f) = [ (λ + 1/4) * fp^4 ]^λ / Γ(λ) * [ Hs^2 / 4 * f^(4λ+1) ] * exp[ -(λ + 1/4) * (fp/f)^4 ]
"""
function compute_ochi_component(f::Real, Hs::Real, fp::Real, λ::Real)
    L = λ + 0.25
    
    # Pre-calculate powers for efficiency
    f4 = f^4
    fp4 = fp^4
    
    # Numerator part
    term1 = (L * fp4)^λ
    # Denominator part (f domain)
    term2 = gamma(λ) * f^(4.0*λ + 1.0)
    # Exponential part
    term3 = exp( -L * (fp4 / f4) )
    
    return (term1 / term2) * (Hs^2 / 4.0) * term3
end