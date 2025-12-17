# TMA.jl

struct TMA{T<:Real} <: AbstractSpectrum
    Hs::T
    Tp::T
    γ::T
    h::T  # Water depth
    Ag::T
end

# TMA Depth Factor (Kitaigorodskii et al. approximation)
function tma_depth_factor(ω, h)
    g = 9.81
    ωh = ω * sqrt(h / g)
    if ωh <= 1.0
        return 0.5 * ωh^2
    elseif ωh < 2.0
        return 1.0 - 0.5 * (2.0 - ωh)^2
    else
        return 1.0
    end
end

function get_density(s::TMA, ω::Real)
    # 1. Get the base JONSWAP density
    # (Usually TMA is normalized differently, but this is the standard approach)
    S_jonswap = get_density(JONSWAP(s.Hs, s.Tp, s.γ, s.Ag), ω)
    
    # 2. Apply depth factor
    return S_jonswap * tma_depth_factor(ω, s.h)
end