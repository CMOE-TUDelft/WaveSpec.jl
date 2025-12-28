#module run_wave


using   DrWatson
@quickactivate "WaveSpec"

using   WaveSpec
using   .Constants
using   .Jonswap
using   .WaveTimeSeries


h0 = 50 #still-water depth
H₀ = 7.11 #m
Tₚ = 12 #s
nω = 257

ω, S, A = jonswap(H₀, Tₚ; 
    plotflag=false, nω = nω)

k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω; seed=100)

sp = SpecStruct( h0, ω, S, A, k, α; Hs = H₀, Tp = Tₚ )


t=0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)
@show η, ϕ, u, w

t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)

t=0.6
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)
@show η, ϕ, u, w
η, px, py = waveAiry1D_pPos(sp, t, 0.1, -0.1)
@show η, px, py

t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)

η, px, py = waveAiry1D_pPos(sp, t, 0.1, -0.1)

return nothing


# end