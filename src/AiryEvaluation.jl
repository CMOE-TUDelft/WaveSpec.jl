using .PhysicalConstants: g

export evaluate_ϕ
export evaluate_u
export evaluate_v
export evaluate_w
export evaluate_η

"""
    _phase(comps::WaveComponents, n::Int, x::Real, y::Real, t::Real)

Compute the phase of the nth Airy wave component at the given spatial
and temporal coordinates.
"""
@inline function _phase(comps::WaveComponents, n::Int, x::Real, y::Real, t::Real)
    comps.kx[n] * x +
    comps.ky[n] * y -
    comps.ω[n] * t +
    comps.phase[n]
end

"""
    evaluate_ϕ(realization::AiryRealization, x::Real, y::Real, z::Real, t::Real)

Evaluate the velocity potential ϕ at a single point and time using a
precomputed Airy realization. The formula for the velocity potential ϕ is given by:

```math
ϕ(x, y, z, t) = \\sum_{i} A_i \\frac{g}{\\omega_i} \\frac{\\cosh(k_i (z + h))}{\\cosh(k_i h)} \\sin(k_{x,i} x + k_{y,i} y - \\omega_i t + \\phi_i)
```

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_ϕ(realization::AiryRealization{T}, x::Real, y::Real, 
                    z::Real, t::Real) where {T}

    comps = realization.components
    h = realization.h

    ϕ = zero(T)
    @inbounds @simd for i in eachindex(comps.ω)
      k = realization.k[i]
      ψ = _phase(comps, i, x, y, t)
      kh = k * h
      coeff =
        if kh < 20.0
            g / comps.ω[i] * cosh(k * (z + h)) / cosh(kh)
        else
            g / comps.ω[i] * exp(k * z)
        end
      ϕ += comps.amplitude[i] * coeff * sin(ψ)
    end
    return ϕ
end

"""
    evaluate_u(realization::AiryRealization, x::Real, y::Real, z::Real, t::Real)

Evaluate the horizontal velocity component u at a single point and time using a
precomputed Airy realization. The formula for the horizontal velocity component u is given by:

```math
u(x, y, z, t) = \\sum_{i} A_i \\omega_i \\frac{k_{x,i}}{k_i} \\frac{\\cosh(k_i (z + h))}{\\sinh(k_i h)} \\cos(k_{x,i} x + k_{y,i} y - \\omega_i t + \\phi_i)
```

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_u(realization::AiryRealization{T}, x::Real, y::Real,
                    z::Real, t::Real) where {T}

    comps = realization.components
    h = realization.h
    u = zero(T)

    @inbounds @simd for i in eachindex(comps.ω)
      k = realization.k[i]
      ψ = _phase(comps, i, x, y, t)
      kh = k * h
      coeff =
        if kh < 20.0
            cosh(k * (z + h)) / sinh(kh)
        else
            exp(k * z)
        end
      cosθ = comps.kx[i] / k
      u += comps.amplitude[i] * comps.ω[i] * cosθ * coeff * cos(ψ)
    end
    return u
end

"""
    evaluate_v(realization::AiryRealization, x::Real, y::Real, z::Real, t::Real)

Evaluate the horizontal velocity component v at a single point and time using a
precomputed Airy realization. The formula for the horizontal velocity component v is given by:

```math
v(x, y, z, t) = \\sum_{i} A_i \\omega_i \\frac{k_{y,i}}{k_i} \\frac{\\cosh(k_i (z + h))}{\\sinh(k_i h)} \\cos(k_{x,i} x + k_{y,i} y - \\omega_i t + \\phi_i)
```

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_v(realization::AiryRealization{T}, x::Real, y::Real,
                    z::Real, t::Real) where {T}

    comps = realization.components
    h = realization.h
    v = zero(T)

    @inbounds @simd for i in eachindex(comps.ω)
      k = realization.k[i]
      ψ = _phase(comps, i, x, y, t)
      kh = k * h
      coeff =
        if kh < 20.0
            cosh(k * (z + h)) / sinh(kh)
        else
            exp(k * z)
        end
      sinθ = comps.ky[i] / k
      v += comps.amplitude[i] * comps.ω[i] * sinθ * coeff * cos(ψ)
    end
    return v
end

"""
    evaluate_w(realization::AiryRealization, x::Real, y::Real, z::Real, t::Real)

Evaluate the vertical velocity component w at a single point and time using a
precomputed Airy realization. The formula for the vertical velocity component w is given by:

```math
w(x, y, z, t) = \\sum_{i} A_i \\omega_i \\frac{\\sinh(k_i (z + h))}{\\sinh(k_i h)} \\sin(k_{x,i} x + k_{y,i} y - \\omega_i t + \\phi_i)
```

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_w(realization::AiryRealization{T}, x::Real, y::Real, 
                    z::Real, t::Real) where {T}

    comps = realization.components
    h = realization.h
    w = zero(T)
    @inbounds @simd for i in eachindex(comps.ω)
      k = realization.k[i]
      ψ = _phase(comps, i, x, y, t)
      kh = k * h
      coeff =
        if kh < 20.0
            sinh(k * (z + h)) / sinh(kh)
        else
            exp(k * z)
        end
      w += comps.amplitude[i] * comps.ω[i] * coeff * sin(ψ)
    end
    return w
end

"""
    evaluate_η(realization, x, y, t)

Evaluate free-surface elevation η at a single point and time using a
precomputed Airy realization. It assumes that the free surface is located
at z = 0. The formula for the free-surface elevation η is given by:

```math
η(x, y, t) = \\sum_{i} A_i \\cos(k_{x,i} x + k_{y,i} y - \\omega_i t + \\phi_i)
```

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_η( realization::AiryRealization{T}, x::Real, y::Real, t::Real) where {T}

    comps = realization.components
    η = zero(T)
    @inbounds @simd for i in eachindex(comps.ω)
        ψ = _phase(comps, i, x, y, t)
        η += comps.amplitude[i] * cos(ψ)
    end
    return η
end
