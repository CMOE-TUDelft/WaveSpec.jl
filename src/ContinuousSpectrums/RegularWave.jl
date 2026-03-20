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
Implements a Dirac delta: returns H^2/8 at f = 1/T, 0 elsewhere.
"""
function get_density(s::RegularWave, f::Real)
    f0 = 1 / s.T
    H = s.H
    if isapprox(f, f0; atol=eps(f0))
        return H^2 / 8  # All energy at f0
    else
        return 0.0
    end
end
