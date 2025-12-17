module WaveFreqDom
# -------------------#
# This next section compiles the wave from 
# all seperate waves into one irregular wave

# using DrWatson
# @quickactivate "WaveSpec"

using WaveSpec.Constants

export waveAiry1DFreq

#----------- waveAiry1DFreq ----------#
function waveAiry1DFreq(h0::Real, ω, A, k, α;
  x::Real=0, z::Real=0, x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1DFreq( params )
  # Ref() is to ensure that broadcast operation is only
  # applied to `t`    

  η, ϕ, u, w = res  

  return η, ϕ, u, w
end


function waveAiry1DFreq( params )    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) + α

  η = A * exp(im*αₘ)

  ϕ = ifelse( ω > 0,
    -im * g* A / ω  * cosh(k*(h0+z)) / cosh(k*h0) * exp(im*αₘ),
    0.0 ) 
  
  u = ifelse( ω > 0, 
    A * ω * cosh(k*(h0+z)) / sinh(k*h0) * exp(im*αₘ),
    0.0 )
  
  w = ifelse( ω > 0, 
    -im * A * ω * sinh.(k*(h0+z)) / sinh(k*h0) * exp(im*αₘ),
    0.0 ) 

  return η, ϕ, u, w
end
#--------- End waveAiry1DFreq --------#


end 