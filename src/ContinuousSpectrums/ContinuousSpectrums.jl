module ContinuousSpectrums

using ..PhysicalConstants
using ..Integration
using ..Truncation
using RecipesBase

export AbstractSpectrum
export JONSWAP, TMA, Bretschneider, Donelan, OchiHubble 
export get_Hs, get_Tp, get_density, get_fmax, get_fmin
export get_phi, compute_ochi_component
export integrate, get_cumulative_energy

# All spectrum models which are defined in this folder/module will inherit from abstract type AbstractSpectrum.
abstract type AbstractSpectrum end

# --- DEFINE AbstractSpectrum's FUNCTIONS INTERFACES ---

"""
    get_Hs(s::AbstractSpectrum)
Returns the significant wave height Hs [m] of the spectrum.
"""
function get_Hs(s::AbstractSpectrum)
    error("Function 'get_Hs' not implemented for $(typeof(s))")
end

"""
    get_Tp(s::AbstractSpectrum)
Returns the peak period Tp [s] of the spectrum.
"""
function get_Tp(s::AbstractSpectrum)
    error("Function 'get_Tp' not implemented for $(typeof(s))")
end

"""
    get_density(s::AbstractSpectrum, f::Real)

Returns the spectral density S(f) in [m²s]. 
Must be implemented by all subtypes.
"""
function get_density(s::AbstractSpectrum, f::Real)
    error("Function 'get_density' not implemented for $(typeof(s))")
end

# --- INCLUDE SPECIFIC SPECTRUM MODELS ---

# Include the specific models
include("JONSWAP.jl")
include("TMA.jl")
include("Bretschneider.jl")
include("Donelan.jl")
include("OchiHubble.jl")


# --- SPECTRAL RANGE ---

"""
    get_fmax(s::AbstractSpectrum; multiplier=5.0)

Returns a recommended maximum frequency for numerical integration 
based on the peak frequency: f_max = multiplier * f_p.
According to literature, f_max between 3 to 5 times f_p is sufficient to capture more than 99% of the energy.
"""
function get_fmax(s::AbstractSpectrum; multiplier=5.0)
    return multiplier / get_Tp(s)
end

function get_fmin(s::AbstractSpectrum; multiplier=1e-4)
    return multiplier / get_Tp(s)
end

# ----------------------

# --- SPECTRAL INTEGRATION ---

""" 
    integrate(s::AbstractSpectrum, fmin, fmax; npoints, method=:symbol, order)

Returns the numeric integral of the spectrum between fmin and fmax, using npoints-1 bins 
and the selected method of integration 
"""
function integrate(s::AbstractSpectrum, fmin::Real, fmax::Real; npoints::Int=100, method::Symbol=:gauss, order::Int=2) 
    if method == :gauss
        return IntegrateGaussQuad(f -> get_density(s, f); order = order, a = fmin, b = fmax, n = npoints)
    elseif method == :trapezoidal 
        return IntegrateTrapezoidal(f -> get_density(s, f); a = fmin, b = fmax, n = npoints)
    else
        throw(ArgumentError("Integration method :$method not recognized. Use :gauss or :trapezoidal."))
    end
end

# Integrate entire spectrum 
function integrate(s::AbstractSpectrum; npoints::Int=200, method::Symbol=:gauss, order::Int=2) 
    if method == :gauss
        return IntegrateGaussQuad(f -> get_density(s, f); order = order, a = get_fmin(s), b = get_fmax(s), n = npoints)
    elseif method == :trapezoidal 
        return IntegrateTrapezoidal(f -> get_density(s, f); a = get_fmin(s), b = get_fmax(s), n = npoints)
    else
        throw(ArgumentError("Integration method :$method not recognized. Use :gauss or :trapezoidal."))
    end
end

"""
    get_cumulative_energy(s::AbstractSpectrum, f::Real, npoints, method=:symbol, order)

Compute the CDF of abstract spectrum at frequency f, using npoints-1 bins 
and the selected method of integration. 
"""
function get_cumulative_energy(s::AbstractSpectrum, f::Real; npoints::Int=100, method::Symbol=:gauss, order::Int=2)
    if method == :gauss
        return IntegrateGaussQuad(f -> get_density(s, f); order = order, a = get_fmin(s), b = f, n = npoints)
    elseif method == :trapezoidal 
        return IntegrateTrapezoidal(f -> get_density(s, f); a = get_fmin(s), b = f, n = npoints)
    else
        throw(ArgumentError("Integration method :$method not recognized. Use :gauss or :trapezoidal."))
    end
end

# ---------------------------

# --- Plots and Visualisation ---
# Plot using light package RecipesBase:
#   1. import Plots in script
#   2. Call this function as plot(DiscreteSpectralSpreading)

@recipe function f(s::AbstractSpectrum; n_points=500)
    # Plot Attributes
    title  --> "Continuous Spectrum PDF"
    xlabel --> "Frequency f [Hz]"
    ylabel --> "Spectral Density S(f) [m²s]"
    grid   --> false

    # User provided features, otherwise default values
    color_spectrum = get(plotattributes, :color, :black)
    spectrum_label = get(plotattributes, :label, "")
    alpha = get(plotattributes, :alpha, 0.3)

    # Get plot range 
    f_range = range(get_fmin(s), get_fmax(s), length=n_points)

    @series begin
        label := spectrum_label
        seriestype := :path
        fillrange := 0
        fillalpha := alpha
        linecolor := color_spectrum
        f_range, [get_density(s, f) for f in f_range]
    end

end

end # module