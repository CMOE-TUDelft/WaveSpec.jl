module Wave
using Random
using Distributions
using Roots
using Plots

export waveparameters

function waveparameters(f, fpeak, d)      # Function that generates the input for a random phase amplitude model, according to the Jonswap spectrum
  g = 9.81
  #CHECK
  E_jon(f,α_la,fpeak,γ,σ) = α_la*g^2*(2*π)^(-4)*f^(-5)*exp(-5/4*(f/fpeak)^(-4))*γ^(exp(-0.5*((f/fpeak)/σ)^2))   #Formulation of the jonswap spectrum

  df   = f[2]-f[1]
  α_la = 0.076*fpeak^0.22                 # According to oceanwiki
  #α_la = 0.0317*fpeak^0.67               # See holthuis pg.162
  γ    = 3.3                              #5.870*fpeak^0.86             # See holthuis pg.162
  σ_a  = 0.07                             #0.0547*fpeak^0.32
  σ_b  = 0.09                             #0.0783*fpeak^0.16

  Random.seed!(1234)
  α =rand(Uniform(0,2π), length(f))        # Generates a random phase difference
  μ = zeros(length(f))
  k = zeros(length(f))
  jon = zeros(length(f))
  m₀ = 0
  m₁ = 0

  for i in 1:length(f)
      if f[i] < fpeak
          σ = σ_a
      else
          σ = σ_b
      end
      #α[i] =  rand(Uniform(0, 2*π))                       # Generation of the random phase
      μ[i] =  √(E_jon(f[i],α_la,fpeak,γ,σ)*df*2)           # Transformation of the Jonswap variance density spectrum to an amplitude spectrum
      jon[i] = E_jon(f[i],α_la,fpeak,γ,σ)
      m₀ += E_jon(f[i],α_la,fpeak,γ,σ)*df
      m₁ += E_jon(f[i],α_la,fpeak,γ,σ)*f[i]*df

      G(K) = √(g*K*tanh(K*d)) - 2π*f[i]     # Formula to compute the wave number (see tutorial from S. van hoof)
      k[i] = abs(find_zero(G, 0.05))        # wave number
  end

  ηₛ = 4*√(m₀)                        # significant wave height
  fₛ = m₁/m₀                          # significant wave frequency

  H(K) = √(g*K*tanh(K*d)) - 2*π*fₛ    # Formula to compute the significant wave number (see tutorial from S. van hoof)
  kₛ = abs(find_zero(H, (0.001, 15), Bisection()))          # significant wave number

  return α, μ, k, kₛ, fₛ, ηₛ, jon
end

###-------------------###
### This next section compiles the wave from all seperate waves into one irregular wave ###

function wave_continous(x, t, f)                      # Function that sums the random phase amplitude equations and is used for the calculation
    Eta = 0
    Phi = 0
    Vxin = 0
    Vzin = 0

    for (i) in 1:length(f)
        Eta  += eta[i]*cos(2*π*f[i]*t - k[i]*x[1]+α[i])
        #CHECK
        Phi  += f[i]*2*π*eta[i]/k[i]* cosh(k[i]*(d-x[2])) /sinh(k[i]*d) *cos(2*pi*f[i]*t - k[i]*x[1] +α[i])
        Vxin += f[i]*2*π*eta[i]     * cosh(k[i]*(d-x[2])) /sinh(k[i]*d) *sin(2*pi*f[i]*t - k[i]*x[1] +α[i])
        Vzin += f[i]*2*π*eta[i]     * sinh(k[i]*(d-x[2])) /sinh(k[i]*d) *cos(2*pi*f[i]*t - k[i]*x[1] +α[i])
    end
    return Eta, Phi, Vxin, Vzin 
end
 
function wave_discrete(x,t,f)                         # Function that discretises the amplitude equation for t and is used for plotting
    Eta = zeros(length(t)) 
    Phi = zeros(length(t))
    Vxin = zeros(length(t))
    Vzin = zeros(length(t))
    for j in 1:length(t)
        Eta[j] = wave_continous(x,t[j],f)[1] 
        Phi[j] = wave_continous(x,t[j],f)[2]
        Vxin[j] = wave_continous(x,t[j],f)[3]
        Vzin[j] = wave_continous(x,t[j],f)[4]
    end
    return Eta, Phi, Vxin, Vzin
end

   # Irregular wave properties
   ηᵢₙ(x,t) = wave_continous(x, t, f)[1]
   ϕᵢₙ(x,t) = wave_continous(x, t, f)[2]
   vxᵢₙ(x,t) = wave_continous(x, t, f)[3]
   vzᵢₙ(x,t) = wave_continous(x, t, f)[4]
   #ϕᵢₙ_∇ₙ_sq(x,t) = (vzᵢₙ(x,t)) .* ϕᵢₙ(x,t)
   #η_sq(x,t) = ηᵢₙ(x,t) .* ηᵢₙ(x,t) 
 
   
   ηᵢₙ(t::Real)   = x -> ηᵢₙ(x,t)
   ϕᵢₙ(t::Real)   = x -> ϕᵢₙ(x,t)
   vxᵢₙ(t::Real)  = x -> vxᵢₙ(x,t)
   vzᵢₙ(t::Real)  = x -> vzᵢₙ(x,t)

end
