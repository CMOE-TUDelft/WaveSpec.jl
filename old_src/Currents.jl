module Currents

# using DrWatson
# @quickactivate "WaveSpec"

using WaveSpec.Constants
using PCHIPInterpolation

export CurrentStat, getCurrent


"""
Struct: CurrentStat
===============

"""
# ---------------------Start---------------------
struct CurrentStat
  h0::Real
  d::Vector{Real} #d=0 at mean sea level
  ux::Vector{Real}
  itp::PCHIPInterpolation.Interpolator

  function CurrentStat(h0::Real, 
    d::Vector{<:Real}, ux::Vector{<:Real})

    # d = append!([-h0], d)
    # ux = append!([0.0], ux)

    # d = append!(d, [0.10*h0])
    # ux = append!(ux, [ux[end]])
    
    # @show d
    # @show ux
    
    # PCHIP Inteprolation
    # itp = Interpolator(d, ux) #no extrapolation
    itp = Interpolator(d, ux, extrapolate = true)
    @show itp
    
    new(h0, d, ux, itp)  
  end

end
# ----------------------End----------------------




"""
Function: getCurrent()
=======================

"""
function getCurrent(cur::CurrentStat, z)

  zp = max(-cur.h0, z)
  zp = min(zp, 0.0)

  return cur.itp(zp)

end
# ---------------------Start---------------------

end