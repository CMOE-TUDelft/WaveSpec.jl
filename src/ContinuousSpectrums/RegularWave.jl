"""
    RegularWave(H, T)

Defines a regular (monochromatic) wave spectrum with:
    - wave height H (m)
    - wave period T (s)
Implements AbstractSpectrum interface.
"""
struct RegularWave{T<:Real} <: AbstractSpectrum
    H::T
    T::T
end

"""
    get_Hs(s::RegularWave)
Returns the wave height (H) for the regular wave.
"""
get_Hs(s::RegularWave) = s.H


"""
    get_Tp(s::RegularWave)
Returns the wave period (T) for the regular wave.
"""
get_Tp(s::RegularWave) = s.T

"""
    get_density(s::RegularWave, f::Real)
Returns the spectral density S(f) for the regular wave.
All energy is at f = 1/T, so S(f) is a narrow Gaussian centered at f0 = 1/T.
"""
function get_density(s::RegularWave, f::Real)
    f0 = 1 / s.T
    H = s.H
    # Total energy = H^2/8, so area under S(f) = H^2/8
    # Use a narrow Gaussian for numerical purposes
    σ = 0.01 * f0  # 1% of central frequency
    norm = (H^2 / 8) / (σ * sqrt(2π))
    return norm * exp(-0.5 * ((f - f0)/σ)^2)
end
