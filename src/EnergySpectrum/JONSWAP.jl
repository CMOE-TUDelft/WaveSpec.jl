# JONSWAP.jl

"""
    JONSWAP(Hs, Tp, γ=nothing)

Creates a JONSWAP spectrum instance with given:
    - significant wave height Hs (m)
    - peak period Tp (s)
    - optional peak enhancement factor γ (if not provided, estimated based on Hs and Tp)
"""


struct JONSWAP{T<:Real} <: AbstractSpectrum
    Hs::T
    Tp::T
    γ::T
    Ag::T
end

function JONSWAP(Hs::T, Tp::T, γ::Union{T, Nothing}=nothing) where {T<:Real}
    actual_γ = isnothing(γ) ? estimate_γ(Hs, Tp) : γ
    Ag = compute_jonswap_norm(actual_γ)
    return JONSWAP{T}(Hs, Tp, actual_γ, Ag)
end

function get_density(s::JONSWAP, ω::Real)
    (ω <= 0) && return 0.0
    ωp = 2π / s.Tp
    ωr = ω / ωp
    σ = ω <= ωp ? 0.07 : 0.09
    
    r = exp(-0.5 * ((ωr - 1) / σ)^2)
    γf = s.γ^r
    
    multiplier = ((s.Hs / 4)^2 / ωp) * s.Ag
    return multiplier * 5.0 * ωr^-5 * exp(-1.25 * ωr^-4) * γf
end


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
    compute_jonswap_norm(γ)
Calculates the factor Ag needed to ensure m₀ = (Hs/4)².
We use a high-resolution quadrature for this pre-computation.
"""
function compute_jonswap_norm(γ::Real)
    # Integration of the non-dimensional shape
    f_shape(wr) = 5.0 * wr^-5 * exp(-1.25 * wr^-4) * γ^exp(-0.5 * ((wr-1)/(wr <= 1 ? 0.07 : 0.09))^2)
    
    # Integrate wr from ~0 to 10 (convergence)
    wr_axis = range(0.01, 10.0, length=1000)
    dwr = wr_axis[2] - wr_axis[1]
    vals = [f_shape(wr) for wr in wr_axis]
    
    area = gaussQuad1D(vals, dwr) # Using your custom module
    return 1.0 / area
end