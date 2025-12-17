module WaveTimeSeries
# -------------------#
# This next section compiles the wave from 
# all seperate waves into one irregular wave

# using DrWatson
# @quickactivate "WaveSpec"

using WaveSpec.Constants
using Parameters


export waveAiry1D
export TimeRampType, timeRamp, timeRampDt
export waveAiry1D_pPos, waveAiry1D_eta
export waveAiry1D_vel


#------------- waveAiry1D ------------#

function waveAiry1D(sp::SpecStruct, t, x, z; 
  x0::Real = 0)    

  return waveAiry1D(sp.h0, sp.ω, sp.A, sp.k, sp.α, 
    t, x, z; 
    x0 = x0)
  
end


function waveAiry1D(h0::Real, ω, A, k, α, 
  t::AbstractArray{<:AbstractFloat},
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D.( Ref(params), t)
  # Ref() is to ensure that broadcast operation is only
  # applied to `t`    

  η = [r for (r, nothing, nothing, nothing) in res]
  ϕ = [r for (nothing, r, nothing, nothing) in res]
  u = [r for (nothing, nothing, r, nothing) in res]
  w = [r for (nothing, nothing, nothing, r) in res]    

  return η, ϕ, u, w
end


function waveAiry1D(h0::Real, ω, A, k, α, t::Real,
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D( params, t)  
  
  η, ϕ, u, w = res

  return η, ϕ, u, w
end


# function waveAiry1D(h0::Real, ω, A, k, α, t,
#   x, z; x0::Real=0)    

#   params = (h0, ω, A, k, α, x0, x, z)
#   res = waveAiry1D.( Ref(params), t)
#   # Ref() is to ensure that broadcast operation is only
#   # applied to `t`    

#   if(res isa Array)
#     η = [r for (r, nothing, nothing, nothing) in res]
#     ϕ = [r for (nothing, r, nothing, nothing) in res]
#     u = [r for (nothing, nothing, r, nothing) in res]
#     w = [r for (nothing, nothing, nothing, r) in res]
  
#   else
#     η, ϕ, u, w = res
#   end

#   return η, ϕ, u, w
# end


function waveAiry1D( params, t::Real)    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) - ω*t + α

  snkh0 = sinh.(k*h0)

  η = sum( A .* cos.(αₘ) )

  ϕ = sum( ifelse.( ω .> 0,
    g* A ./ ω .* sin.(αₘ) .* cosh.(k*(h0+z)) ./ cosh.(k*h0),
    0.0 ) )
  
  u = sum( ifelse.( ω .> 0, 
    A .* ω .* cos.(αₘ) .* cosh.(k*(h0+z)) ./ snkh0,
    0.0 ) )
  
  w = sum( ifelse.( ω .> 0, 
    A .* ω .* sin.(αₘ) .* sinh.(k*(h0+z)) ./ snkh0,
    0.0 ) )

  return η, ϕ, u, w
end
#----------- End waveAiry1D ----------#



#------------- waveAiry1D_eta ------------#

function waveAiry1D_eta(sp::SpecStruct, t, x, z; 
  x0::Real = 0)    

  return waveAiry1D_eta(sp.h0, sp.ω, sp.A, sp.k, sp.α, 
    t, x, z; 
    x0 = x0)
  
end


function waveAiry1D_eta(h0::Real, ω, A, k, α, t::Real,
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  
  return waveAiry1D_eta( params, t)  
  
end


function waveAiry1D_eta( params, t::Real)    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) - ω*t + α

  η = sum( A .* cos.(αₘ) )

  return η
end
#----------- End waveAiry1D_eta ----------#



#------------- waveAiry1D_pPos ------------#

# function waveAiry1D_pPos(sp::SpecStruct, t, x, z; 
#   x0::Real = 0)    

#   return waveAiry1D_pPos(sp.h0, sp.ω, sp.A, sp.k, sp.α, 
#     t, x, z; 
#     x0 = x0)
  
# end
function waveAiry1D_pPos(sp::SpecStruct, 
  t::Real, x::Real, z::Real; 
  x0::Real = 0)    

  @unpack h0, k, ω, ηVec, pxVec, pzVec = sp

  αm = k*(x-x0) - ω*t 
  expiαm = exp.(im*αm)
  coefz = k*(h0+z)

  η = real( transpose(ηVec) * expiαm )
  px = transpose(cosh.(coefz)) * real( (pxVec .* expiαm) )
  pz = transpose(sinh.(coefz)) * real( (pzVec .* expiαm) )
  
  return η, px, pz
end


function waveAiry1D_pPos(sp::SpecStruct, 
  t::AbstractArray{<:AbstractFloat}, 
  x::Real, z::Real; 
  x0::Real = 0)    

  res = waveAiry1D_pPos.(Ref(sp), t, Ref(x), Ref(z); x0 = x0)  

  η = [r for (r, nothing, nothing) in res]
  px = [r for (nothing, r, nothing) in res]
  py = [r for (nothing, nothing, r) in res]  

  return η, px, py
end


function waveAiry1D_pPos(h0::Real, ω, A, k, α, 
  t::AbstractArray{<:AbstractFloat},
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D_pPos.( Ref(params), t)
  # Ref() is to ensure that broadcast operation is only
  # applied to `t`    

  η = [r for (r, nothing, nothing) in res]
  px = [r for (nothing, r, nothing) in res]
  py = [r for (nothing, nothing, r) in res]  

  return η, px, py
end


function waveAiry1D_pPos(h0::Real, ω, A, k, α, t::Real,
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D_pPos( params, t)  
  
  η, px, py = res

  return η, px, py
end


function waveAiry1D_pPos( params, t::Real)    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) - ω*t + α

  snkh0 = sinh.(k*h0)

  η = sum( A .* cos.(αₘ) )  
  
  px = sum( ifelse.( ω .> 0, 
    - A .* sin.(αₘ) .* cosh.(k*(h0+z)) ./ snkh0,
    0.0 ) )
  
  py = sum( ifelse.( ω .> 0, 
    A .* cos.(αₘ) .* sinh.(k*(h0+z)) ./ snkh0,
    0.0 ) )

  return η, px, py
end
#----------- End waveAiry1D_pPos ----------#



#------------- waveAiry1D_vel ------------#

# function waveAiry1D_vel(sp::SpecStruct, t, x, z; 
#   x0::Real = 0)    

#   return waveAiry1D_vel(sp.h0, sp.ω, sp.A, sp.k, sp.α, 
#     t, x, z; 
#     x0 = x0)
  
# end
# function waveAiry1D_vel(sp::SpecStruct, t::Real, x::Real, z::Real; 
#   x0::Real = 0)    

#   @unpack h0, k, ω, uVec, wVec, α, A = sp

#   αm = k*(x-x0) - ω*t 
#   expiαm = exp.(im*αm)
#   expiα = exp.(im*α)
#   coefz = k*(h0+z)
#   ph = expiαm .* expiα
#   cosαm = real.(ph)
#   sinαm = imag.(ph)

#   u = (uVec .* cosh.(coefz))' * cosαm
#   w = (wVec .* sinh.(coefz))' * sinαm

#   return u, w    
# end
function waveAiry1D_vel(sp::SpecStruct, t::Real, x::Real, z::Real; 
  x0::Real = 0)    

  @unpack h0, k, ω, uVec, wVec = sp

  αm = k*(x-x0) - ω*t 
  expiαm = exp.(im*αm)
  coefz = k*(h0+z)

  # u = real(transpose(cosh.(coefz)) * (uVec .* expiαm) )
  # w = real(transpose(sinh.(coefz)) * (wVec .* expiαm) ) 

  u = transpose(cosh.(coefz)) * real( (uVec .* expiαm) )
  w = transpose(sinh.(coefz)) * real( (wVec .* expiαm) ) 

  return u, w    
end


function waveAiry1D_vel(h0::Real, ω, A, k, α, t::Real,
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  return waveAiry1D_vel( params, t)  
  
end


function waveAiry1D_vel( params, t::Real)    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) - ω*t + α

  snkh0 = sinh.(k*h0)

  u = sum( ifelse.( ω .> 0, 
    A .* ω .* cos.(αₘ) .* cosh.(k*(h0+z)) ./ snkh0,
    0.0 ) )
  
  w = sum( ifelse.( ω .> 0, 
    A .* ω .* sin.(αₘ) .* sinh.(k*(h0+z)) ./ snkh0,
    0.0 ) )

  return u, w
end
#----------- End waveAiry1D_vel ----------#




#----------- TimeRampType ------------#
abstract type TimeRampTypeAll end
struct TimeRampType1 <: TimeRampTypeAll
  t0::Real
  t1::Real  
end

struct TimeRampType2 <: TimeRampTypeAll
  t0::Real
  t1::Real  
  t2::Real  
  tEnd::Real  
end

function TimeRampType(t0, t1)
  TimeRampType1(t0, t1)
end

function TimeRampType(t0, t1, t2, tEnd)
  TimeRampType2(t0, t1, t2, tEnd)
end

function timeRamp(t::AbstractArray{<:AbstractFloat}, 
  tRampObj::TimeRampTypeAll)
  
  timeRamp.(t, Ref(tRampObj))
end

function timeRamp(t::Real, tRampObj::TimeRampType1)
  timeRamp(t, tRampObj.t0, tRampObj.t1)
end

function timeRamp(t::Real, tRampObj::TimeRampType2)
  timeRamp(t, 
    tRampObj.t0, tRampObj.t1, 
    tRampObj.t2, tRampObj.tEnd)
end

function timeRampDt(t::Real, tRampObj::TimeRampType1)
  timeRampDt(t, tRampObj.t0, tRampObj.t1)
end
#--------- End TimeRampType ----------#




#------------- timeRamp --------------#
function timeRamp(t::Real, t0, t1, t2, tEnd)

  if(t1 ≤ t ≤ t2)
    return 1.0
  elseif(t0 < t < t1)
    return 0.5*( 1.0 - cos(π*(t - t0) / (t1-t0)) )
  elseif(t2 < t < tEnd)
    return 0.5*( 1.0 - cos(π*(tEnd - t) / (tEnd-t2)) )
  else
    return 0.0
  end

end

function timeRamp(t::Real, t0, t1)

  if(t1 ≤ t )
    return 1.0
  elseif(t0 < t < t1)
    return 0.5*( 1.0 - cos(π*(t - t0) / (t1-t0)) )  
  else
    return 0.0
  end

end


function timeRampDt(t::Real, t0, t1)

  if(t1 ≤ t )
    return 0.0
  elseif(t0 < t < t1)
    return 0.5*π/(t1-t0) * sin(π*(t - t0) / (t1-t0))    
  else
    return 0.0
  end

end
#----------- End timeRamp ------------#


end 