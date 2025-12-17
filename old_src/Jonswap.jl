module Jonswap

# using DrWatson
# @quickactivate "WaveSpec"

using Plots
using WaveSpec.Constants

export jonswap


"""
Function: jonswap()
===================

# Description
# ωc = cutoff freq (default 33/Tₚ)
"""
# ---------------------Start---------------------

function jonswap(Hₛ::Real, Tₚ::Real, 
  γ::Real = getJonswapPeakedness(Hₛ, Tₚ);
  τa::Real = 0.07, τb::Real = 0.09,
  plotflag::Bool = false, 
  ωc::Real = 33.0/Tₚ, 
  nω::Int64 = 257, 
  Ag::Real = 0,
  plotloc = "./")


  # Only coded for the case N=5, M=4 
  # N/M = 1.25 (see WAFO jonswap formulation)

  # # Make nω odd for Simpson's integration
  # nω = iseven(nω) ? nω+1 : nω
  
  ωₚ = 2*π/Tₚ
  ω = range(0, ωc, nω)
  dω = ω[2] - ω[1]
  ωᵣ = ω / ωₚ

  # Shape parameter
  # τ = [x <= 1 ? τa : τb for x in ωᵣ]
  τ(ωᵣ) = ifelse.(ωᵣ .<= 1, τa, τb)  

  # Peak amplification factor
  γf = γ .^ exp.( -0.5*( (ωᵣ.-1)./τ(ωᵣ) ).^2 )

  # Calculation of Ag
  if(Ag ≈ 0)
    S2 = ifelse.(ω .> 0 ,
      5 * ωᵣ.^(-5) .* exp.( -1.25 * (ωᵣ).^(-4) ) .* γf,
      0) 

    mask1 = findall(ωᵣ .<= 1)
    mask2 = (mask1[end]:length(ωᵣ))

    sS2 = gaussQuad1D(S2[mask1], dω/ωₚ) + 
      gaussQuad1D(S2[mask2], dω/ωₚ)    
    Ag = 1/sS2
  end

  println("Hₛ \t ",round(Hₛ; digits=3))
  println("Tₚ \t ",round(Tₚ; digits=3))
  println("ωₚ \t ",round(ωₚ; digits=3))
  println("γ \t ",round(γ; digits=3))
  println("Ag \t ",round(Ag; digits=3))
    
  # S = ifelse.(ω .> 0 ,
  #   α * g^2 * ω.^(-5) .* γf .* exp.( -1.25 * (ωᵣ).^(-4) ),
  #   0) 
  
  S = ifelse.(ω .> 0 ,
    ((Hₛ/4)^2/ωₚ.*Ag) * 5 * ωᵣ.^(-5) .* exp.(-1.25*(ωᵣ).^(-4)) .* γf,
    0) 
  
  A = .√( S * dω * 2 )

  # m₀ = simpsonInteg1D(S, dω) 
  m₀ = gaussQuad1D(S, dω) 
  m₁ = gaussQuad1D(ω.*S, dω) 
  cHₛ = 4*√(m₀)
  cωₛ = m₁ / m₀
  println("Computed Hₛ \t ", round(cHₛ; digits=3))
  println("Computed ωₛ \t ", round(cωₛ; digits=3))
  println("Max(S) \t", round(maximum(S); digits=3))
  println("Max(A) \t", round(maximum(A); digits=3))
  

  if(plotflag)
    plt1 = plot(ω, S, linewidth=3, dpi=330, 
    xlabel = "ω (rad/s)",
    ylabel = "S (m^2 s/rad)",
    label = "Spectral density")#,
    #xticks = 0:0.2:ωc,
    #yticks = 0:2:maximum(S))
    #plot!(size = (1280, 640))
    savefig(plt1,plotloc*"jonswap_specdensity.png")

    plt2 = plot(ω, A, linewidth=3, dpi=330,
    xlabel = "ω (rad/s)",
    ylabel = "Amplitude (m)",
    label = "Amplitude")#,
    #xticks = 0:0.2:ωc)
    #plot!(size = (1280, 640))
    savefig(plt2,plotloc*"jonswap_amplitude.png")
  end

  return [ω;], S, A, cHₛ, cωₛ
end


function getJonswapPeakedness(Hₛ::Real, Tₚ::Real)
  
  x = Tₚ / √(Hₛ)
  if(x<5.14285714285714)
    D = 0.036-0.0056*x; 
    # approx 5.061*Hm0^2/Tp^4*(1-0.287*log(gam));
    γ = minimum([exp(3.484*( 1-0.1975*D*x^4 )), 7.0])
  else
    println("ERR: Check criterion Tₚ / √(Hₛ) = ", x)
    γ = 1.0
  end

  return maximum([γ, 1.0])
end


"""
Estimation of Hs, given Tp and gamma
"""

function estimateHs(Tp, γ = 3.3)
  return (0.11661 + 0.01581*γ - 0.00065*γ*γ)*Tp*Tp
end

# ----------------------End----------------------

end