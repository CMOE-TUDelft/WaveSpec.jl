export evaluate_ϕ
export evaluate_u
export evaluate_v
export evaluate_w
export evaluate_η
export evaluate_fields

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
    _vertical_profiles(k, e2kh, z)

Return `(C, S)` with `C = e^{kz} + e^{-k(z+2h)}` and `S = e^{kz} - e^{-k(z+2h)}`,
which are `2e^{-kh} cosh(k(z+h))` and `2e^{-kh} sinh(k(z+h))`. Divided by the
per-component `1 ± e^{-2kh}` factors stored in `AiryRealization`, they give the
exact hyperbolic depth ratios without overflowing for large `kh`, so no
shallow/deep-water branch is needed.

`e^{-k(z+2h)}` is computed as `e2kh / e^{kz}`, so only one `exp` is needed per
component. When both underflow to zero (very deep water near the bed) the term
is set to zero instead of `0/0`.
"""
@inline function _vertical_profiles(k, e2kh, z)
    e⁺ = exp(k * z)
    e⁻ = ifelse(iszero(e⁺), zero(e⁺), e2kh / e⁺)
    return e⁺ + e⁻, e⁺ - e⁻
end

# Accumulator type: promotes the realization precision with the coordinates.
@inline _acctype(::Type{T}, coords::Real...) where {T} = float(promote_type(T, map(typeof, coords)...))

# Threaded point-wise driver shared by all in-place evaluators. `eachindex`
# checks that all arrays share the same axes, which makes `@inbounds` safe.
function _pointwise!(kernel::F, out::AbstractVector, realization::AiryRealization,
                     coords::AbstractVector...) where {F}
    idx = eachindex(out, coords...)
    Threads.@threads for j in idx
        @inbounds out[j] = kernel(realization, map(c -> (@inbounds c[j]), coords)...)
    end
    return out
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
    ϕ = zero(_acctype(T, x, y, z, t))
    @inbounds @simd for i in eachindex(comps.ω)
        C, _ = _vertical_profiles(realization.k[i], realization.e2kh[i], z)
        ϕ += comps.amplitude[i] * realization.cϕ[i] * C * sin(_phase(comps, i, x, y, t))
    end
    return ϕ
end

evaluate_ϕ!(ϕ::AbstractVector, realization::AiryRealization,
            x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
            z::AbstractVector{<:Real}, t::AbstractVector{<:Real}) =
    _pointwise!(evaluate_ϕ, ϕ, realization, x, y, z, t)

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
    u = zero(_acctype(T, x, y, z, t))
    @inbounds @simd for i in eachindex(comps.ω)
        C, _ = _vertical_profiles(realization.k[i], realization.e2kh[i], z)
        u += comps.amplitude[i] * realization.cu[i] * C * cos(_phase(comps, i, x, y, t))
    end
    return u
end

evaluate_u!(u::AbstractVector, realization::AiryRealization,
            x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
            z::AbstractVector{<:Real}, t::AbstractVector{<:Real}) =
    _pointwise!(evaluate_u, u, realization, x, y, z, t)

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
    v = zero(_acctype(T, x, y, z, t))
    @inbounds @simd for i in eachindex(comps.ω)
        C, _ = _vertical_profiles(realization.k[i], realization.e2kh[i], z)
        v += comps.amplitude[i] * realization.cv[i] * C * cos(_phase(comps, i, x, y, t))
    end
    return v
end

evaluate_v!(v::AbstractVector, realization::AiryRealization,
            x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
            z::AbstractVector{<:Real}, t::AbstractVector{<:Real}) =
    _pointwise!(evaluate_v, v, realization, x, y, z, t)

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
    w = zero(_acctype(T, x, y, z, t))
    @inbounds @simd for i in eachindex(comps.ω)
        _, S = _vertical_profiles(realization.k[i], realization.e2kh[i], z)
        w += comps.amplitude[i] * realization.cw[i] * S * sin(_phase(comps, i, x, y, t))
    end
    return w
end

evaluate_w!(w::AbstractVector, realization::AiryRealization,
            x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
            z::AbstractVector{<:Real}, t::AbstractVector{<:Real}) =
    _pointwise!(evaluate_w, w, realization, x, y, z, t)

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
function evaluate_η(realization::AiryRealization{T}, x::Real, y::Real, t::Real) where {T}
    components = realization.components
    η = zero(_acctype(T, x, y, t))
    @inbounds @simd for i in eachindex(components.ω)
        η += components.amplitude[i] * cos(_phase(components, i, x, y, t))
    end
    return η
end

evaluate_η!(η::AbstractVector, realization::AiryRealization,
            x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
            t::AbstractVector{<:Real}) =
    _pointwise!(evaluate_η, η, realization, x, y, t)

"""
    evaluate_fields(realization::AiryRealization, x::Real, y::Real, z::Real, t::Real)

Evaluate all Airy fields at a single point and time in one pass over the
components and return them as a `NamedTuple` `(η, ϕ, u, v, w)`. `η` is the
free-surface elevation (independent of `z`); `ϕ, u, v, w` are evaluated at `z`.

Each component's phase, `sincos` and vertical profile are computed once and
shared by the five fields, so this costs about as much as a single
`evaluate_ϕ` call instead of five separate evaluations.
"""
function evaluate_fields(realization::AiryRealization{T}, x::Real, y::Real,
                         z::Real, t::Real) where {T}
    comps = realization.components
    R = _acctype(T, x, y, z, t)
    η = ϕ = u = v = w = zero(R)
    @inbounds @simd for i in eachindex(comps.ω)
        s, c = sincos(_phase(comps, i, x, y, t))
        C, S = _vertical_profiles(realization.k[i], realization.e2kh[i], z)
        a = comps.amplitude[i]
        aCc = a * C * c
        η += a * c
        ϕ += a * realization.cϕ[i] * C * s
        u += realization.cu[i] * aCc
        v += realization.cv[i] * aCc
        w += a * realization.cw[i] * S * s
    end
    return (η = η, ϕ = ϕ, u = u, v = v, w = w)
end

"""
    evaluate_fields!(η, ϕ, u, v, w, realization, x, y, z, t)

In-place, multi-threaded version of [`evaluate_fields`](@ref) over the points
`(x[j], y[j], z[j], t[j])`. All arrays must share the same axes.
"""
function evaluate_fields!(η::AbstractVector, ϕ::AbstractVector, u::AbstractVector,
                          v::AbstractVector, w::AbstractVector, realization::AiryRealization,
                          x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                          z::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    idx = eachindex(η, ϕ, u, v, w, x, y, z, t)
    Threads.@threads for j in idx
        @inbounds begin
            f = evaluate_fields(realization, x[j], y[j], z[j], t[j])
            η[j] = f.η; ϕ[j] = f.ϕ; u[j] = f.u; v[j] = f.v; w[j] = f.w
        end
    end
    return (η = η, ϕ = ϕ, u = u, v = v, w = w)
end

# --- Gridded evaluation -------------------------------------------------------

# Method of `AiryWaves.generate_sea` (documented there). It lives here because it
# needs `realize` and the kernels above, which are defined after `AiryWaves`.
function AiryWaves.generate_sea(state::AiryWaves.AiryState, x::AbstractArray{<:Real},
                                y::AbstractArray{<:Real}, z::AbstractArray{<:Real},
                                t::AbstractArray{<:Real}; vars = [:η, :ϕ, :u, :v, :w])
    requested = vars isa Symbol ? (vars,) : Tuple(vars)
    if isempty(requested)
        throw(ArgumentError("No valid variables requested. Choose from :η, :ϕ, :u, :v, :w"))
    end
    realization = realize(state)
    res = _generate_sea(realization, vec(x), vec(y), vec(z), vec(t), requested)
    return (; (v => res[v] for v in unique(requested) if haskey(res, v))...)
end

function _generate_sea(realization::AiryRealization{T}, x::AbstractVector, y::AbstractVector,
                       z::AbstractVector, t::AbstractVector, requested) where {T}
    nx, ny, nz, nt = length(x), length(y), length(z), length(t)
    res = Dict{Symbol, Array{T, 4}}()

    if :η in requested
        η = Array{T, 4}(undef, nx, ny, 1, nt)
        Threads.@threads for I in CartesianIndices(η)
            ix, iy, _, it = Tuple(I)
            @inbounds η[I] = evaluate_η(realization, x[ix], y[iy], t[it])
        end
        res[:η] = η
    end

    kinematic = (:ϕ, :u, :v, :w)
    wanted = map(in(requested), kinematic)
    if any(wanted)
        # Unrequested fields get an empty placeholder so the loop stays type-stable.
        ϕ, u, v, w = map(b -> b ? Array{T, 4}(undef, nx, ny, nz, nt) : Array{T, 4}(undef, 0, 0, 0, 0), wanted)
        Threads.@threads for I in CartesianIndices((nx, ny, nz, nt))
            ix, iy, iz, it = Tuple(I)
            @inbounds begin
                f = evaluate_fields(realization, x[ix], y[iy], z[iz], t[it])
                wanted[1] && (ϕ[I] = f.ϕ)
                wanted[2] && (u[I] = f.u)
                wanted[3] && (v[I] = f.v)
                wanted[4] && (w[I] = f.w)
            end
        end
        for (name, arr, b) in zip(kinematic, (ϕ, u, v, w), wanted)
            b && (res[name] = arr)
        end
    end
    return res
end
