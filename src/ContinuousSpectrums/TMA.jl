# TMA.jl

"""
    TMA(Hs, Tp, h, γ=nothing)

The TMA spectrum (Texel-Marsen-Arsloe) is a shallow-water extension 
of the JONSWAP spectrum (mathematically defined as a JONSWAP spectrum 
multiplied by the Kitaigorodskii factor Φ(f, h) ). The Φ function is 
essentially a depth-scaling law that describes how the equilibrium range 
of the spectrum (the high-frequency tail) is limited by water depth. In this
model we will be implementing a Kitaigorodskii factor (Kitaigorodskii 
et al. (1975)) approximation.
"""
struct TMA{T<:Real} <: AbstractSpectrum
    js::JONSWAP{T}  # The underlying JONSWAP model
    h::T            # Water depth [m]
    # We store a specific Ag for TMA because the depth factor changes the area
    Ag_tma::T       
end

# --- 1. Primary Constructor ---
function TMA(Hs::Real, Tp::Real, h::Real, γ::Union{Real, Nothing}=nothing)
    # 1. Create the base JONSWAP (handles γ estimation internally)
    base_js = isnothing(γ) ? JONSWAP(Hs, Tp) : JONSWAP(Hs, Tp, γ)
    
    T = promote_type(typeof(Hs), typeof(Tp), typeof(h))
    
    # 2. Re-calculate normalization because phi(f,h) removes energy
    # We need a new Ag to ensure m0 = (Hs/4)²
    Ag_tma = compute_tma_normalization(base_js, h)
    
    return TMA{T}(base_js, h, Ag_tma)
end

# --- 2. The Phi Factor (Kitaigorodskii) ---
function get_phi(f::Real, h::Real)
    # Dimensionless depth parameter
    ωh = 2π * f * sqrt(h / g)
    
    if ωh <= 1.0
        return 0.5 * ωh^2
    elseif ωh < 2.0
        return 1.0 - 0.5 * (2.0 - ωh)^2
    else
        return 1.0
    end
end

# --- 3. Density Function (The Exploit) ---
function get_density(s::TMA, f::Real)
    # We call the underlying JONSWAP density
    # But we must swap the JONSWAP Ag for our TMA-specific Ag
    js_raw = get_density(s.js, f)
    
    # Since get_density(js) already includes s.js.Ag, 
    # we normalize it out and apply the TMA Ag + Phi factor
    return (js_raw / s.js.Ag) * s.Ag_tma * get_phi(f, s.h)
end

# --- 4. Normalization ---
function compute_tma_normalization(js::JONSWAP, h::Real)
    # Integrate the shape: [JONSWAP_shape] * [Phi]
    # We use js.fp and js.γ from the existing struct
    function tma_shape(f::Real)
        fr = f / js.fp
        # Reuse JONSWAP's internal shape logic
        σ = (fr <= 1.0) ? js.σ₁ : js.σ₂
        b = exp(-0.5 * ((fr - 1.0) / σ)^2)
        S_js_unit = (5.0 * fr^-5 * exp(-1.25 * fr^-4)) * js.γ^b
        
        return (1.0 / js.fp) * S_js_unit * get_phi(f, h)
    end

    total_area = IntegrateGaussQuad(tma_shape, order=2, a=1e-4, b=5.0, n=250)
    return 1.0 / total_area
end