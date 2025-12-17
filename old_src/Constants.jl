module Constants

using Roots: find_zero
using Random
using Distributions

export g, dispersionRel, dispersionRelAng
export randomPhase
export simpsonInteg1D, trapzInteg1D, gaussQuad1D
export SpecStruct


"""
Common
======

"""
# ---------------------Start---------------------

g = 9.81 #accel due to gravity


function randomPhase(ω; seed::Int64 = -1)

  # Seeding is needed to generate the same 
  # array of random numbers every time
  if(seed > 0)
    Random.seed!(seed) 
  end
  ϕ = rand(Uniform(-π,π), length(ω))
  return ϕ
end
# ----------------------End----------------------




"""
Struct: SpecStruct
===============

"""
# ---------------------Start---------------------
struct SpecStruct
  h0::Real
  ω::Vector{Real}
  S::Vector{Real}
  A::Vector{Real}
  k::Vector{Real}
  α::Vector{Real}
  nω::Integer
  Hs::Real 
  Tp::Real  

  ηVec::Vector{Complex}
  ϕVec::Vector{Complex}
  uVec::Vector{Complex}
  wVec::Vector{Complex}
  pxVec::Vector{Complex}  
  pzVec::Vector{Complex}  


  function SpecStruct( h0::Real, ω::Vector{<:Real},
    S::Vector{<:Real}, A::Vector{<:Real}, k::Vector{<:Real},
    α::Vector{<:Real};
    Hs = -99.0, Tp = -99.0 )

    nω = length(ω)    
    ηVec, ϕVec, uVec, wVec, pxVec, pzVec = calcSpecConstants(h0, ω, A, k, α)    
    new( h0, ω, S, A, k, α, nω, Hs, Tp,
      ηVec, ϕVec, uVec, wVec, pxVec, pzVec )    
  end

end


function calcSpecConstants(h0, ω, A, k, α)

  snkh0 = sinh.(k*h0)
  expiα = exp.(im*α)

  η = ifelse.( ω .> 0,
    A .* expiα,
    A .+ 0.0*im )  # To allow for mean value

  ϕ = ifelse.( ω .> 0,
    -im * g* A ./ ω .* expiα ./ cosh.(k*h0),
    0.0 + 0.0*im )

  u = ifelse.( ω .> 0, 
    A .* ω .* expiα ./ snkh0,
    0.0 + 0.0*im )

  w = ifelse.( ω .> 0, 
    -im * A .* ω .* expiα ./ snkh0,
    0.0 + 0.0*im ) 
  
  px = ifelse.( ω .> 0, 
    im * A .* expiα ./ snkh0,
    0.0 + 0.0*im )
  
  pz = ifelse.( ω .> 0, 
    A .* expiα ./ snkh0,
    0.0 + 0.0*im )  
  
  return η, ϕ, u, w, px, pz
end


# function calcSpecConstants(h0, ω, A, k, α)

#   snkh0 = sinh.(k*h0)  

#   η = ifelse.( ω .> 0,
#     A,
#     0.0 )

#   ϕ = ifelse.( ω .> 0,
#     g* A ./ ω ./ cosh.(k*h0),
#     0.0 )

#   u = ifelse.( ω .> 0, 
#     A .* ω ./ snkh0,
#     0.0 + 0.0*im )

#   w = ifelse.( ω .> 0, 
#     A .* ω ./ snkh0,
#     0.0 + 0.0*im ) 
  
#   px = ifelse.( ω .> 0, 
#     A ./ snkh0,
#     0.0 )
  
#   pz = ifelse.( ω .> 0, 
#     A  ./ snkh0,
#     0.0 )

#   return η, ϕ, u, w, px, pz
# end
# ----------------------End----------------------




"""
Function: dispersionRel()
=======================

"""
# ---------------------Start---------------------
function dispersionRel(h::Real, T::Real; msg::Bool=true)

  f(L) = L-(g/2/π*T*T*tanh(2*π/L*h))
  L = find_zero(f,0)

  if(msg)
    println("Time Period \t T \t ",T)
    println("Water Depth \t h \t ",h)
    println("Wave Length \t L \t ",L)
    println("Wave Celerity \t C \t ",L/T)
    println("Disp. Regime \t h/L \t ",round.(h/L; digits=3))
  end
  
  return L
end


function dispersionRelAng(h::Real, ω::Real; msg::Bool=true)

  f(k) = ω^2-(g*k*tanh(k*h))
  k = find_zero(f,0)

  if(msg)
    println("Angular Freq \t ω \t ",ω)
    println("Water Depth \t h \t ",h)
    println("Wave number \t k \t ",k)
    println("Wave Celerity \t C \t ",ω/k)
    println("Disp. Regime \t h/L \t ",round.(h*k/2π; digits=3))
  end
  
  return k
end
# ----------------------End----------------------




"""
Integration
===========

"""
# ---------------------Start---------------------
function simpsonInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 4, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/3.0 * w .* y)
end


function trapzInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 2, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/2.0 * w .* y)
end


function gaussQuad1D(y, dx)
  # Gauss-Legendre Quadrature 2 point
  q1dx = 0.5*(1 - 1/√(3))
  q2dx = 0.5*(1 + 1/√(3))

  ly = y[1:end-1]
  ry = y[2:end]

  yq1 = (ry - ly)*q1dx + ly
  yq2 = (ry - ly)*q2dx + ly

  return sum(yq1 + yq2)*dx/2
end
# ----------------------End----------------------


end