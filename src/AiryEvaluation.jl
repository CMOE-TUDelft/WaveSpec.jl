using .PhysicalConstants: g

export evaluate_phi
export evaluate_u
export evaluate_v
export evaluate_w
export evaluate_eta

@inline function _phase(comps, n, x, y, t)
    comps.kx[n] * x +
    comps.ky[n] * y -
    comps.ω[n] * t +
    comps.phase[n]
end

function evaluate_phi(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    z::Real,
    t::Real,
) where {T}

    comps = realization.components

    h = realization.h

    ϕ = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        k = realization.k[n]

        ψ = _phase(comps, n, x, y, t)

        kh = k * h

        coeff =
            if kh < 20.0
                g / comps.ω[n] *
                cosh(k * (z + h)) /
                cosh(kh)
            else
                g / comps.ω[n] *
                exp(k * z)
            end

        ϕ += comps.amplitude[n] * coeff * sin(ψ)

    end

    return ϕ
end

function evaluate_u(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    z::Real,
    t::Real,
) where {T}

    comps = realization.components

    h = realization.h

    u = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        k = realization.k[n]

        ψ = _phase(comps, n, x, y, t)

        kh = k * h

        coeff =
            if kh < 20.0
                cosh(k * (z + h)) /
                sinh(kh)
            else
                exp(k * z)
            end

        cosθ = comps.kx[n] / k

        u +=
            comps.amplitude[n] *
            comps.ω[n] *
            cosθ *
            coeff *
            cos(ψ)

    end

    return u
end

function evaluate_v(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    z::Real,
    t::Real,
) where {T}

    comps = realization.components

    h = realization.h

    v = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        k = realization.k[n]

        ψ = _phase(comps, n, x, y, t)

        kh = k * h

        coeff =
            if kh < 20.0
                cosh(k * (z + h)) /
                sinh(kh)
            else
                exp(k * z)
            end

        sinθ = comps.ky[n] / k

        v +=
            comps.amplitude[n] *
            comps.ω[n] *
            sinθ *
            coeff *
            cos(ψ)

    end

    return v
end

function evaluate_w(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    z::Real,
    t::Real,
) where {T}

    comps = realization.components

    h = realization.h

    w = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        k = realization.k[n]

        ψ = _phase(comps, n, x, y, t)

        kh = k * h

        coeff =
            if kh < 20.0
                sinh(k * (z + h)) /
                sinh(kh)
            else
                exp(k * z)
            end

        w +=
            comps.amplitude[n] *
            comps.ω[n] *
            coeff *
            sin(ψ)

    end

    return w
end

"""
    evaluate_eta(realization, x, y, t)

Evaluate free-surface elevation η at a single point and time using a
precomputed Airy realization.

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_eta(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    t::Real,
) where {T}

    comps = realization.components

    η = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        ψ =
            comps.kx[n] * x +
            comps.ky[n] * y -
            comps.ω[n] * t +
            comps.phase[n]

        η += comps.amplitude[n] * cos(ψ)

    end

    return η
end
