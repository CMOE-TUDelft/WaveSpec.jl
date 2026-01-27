#module run_wave


using   DrWatson
@quickactivate "WaveSpec"

using   TickTock
using   WaveSpec
using   Plots
using   .Constants
using   .Jonswap
using   .WaveTimeSeries

pltFolder = datadir("output","res_")

h0 = 50 #still-water depth
H₀ = 7.11 #m
Tₚ = 12 #s
nω = 257

ω, S, A = jonswap(H₀, Tₚ; 
    plotflag=true, plotloc=pltFolder, nω = nω)

k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω; seed=100)

sp = SpecStruct( h0, ω, S, A, k, α; Hs = H₀, Tp = Tₚ )

scatter(ω, k)
savefig(pltFolder*"dispersion.png")

scatter(ω, α)
savefig(pltFolder*"phase.png")

t=0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)
@show η, ϕ, u, w

t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)

t=0.6
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)
@show η, ϕ, u, w
η, px, pz = waveAiry1D_pPos(sp, t, 0.1, -0.1)
@show η, px, pz

tick()
@show u, w = waveAiry1D_vel(sp, 10, 0.1, -0.1)
@show u, w = waveAiry1D_vel(sp, 100, 10, -1)
tock()

tick()
for i = 1:1000
    local u, w    
    waveAiry1D_vel(sp, 10, 0.1, -0.1)
    waveAiry1D_vel(sp, 100, 10, -1)
end
tock()


tick()
@show waveAiry1D_eta(sp, 10, 0.1, -0.1)
@show waveAiry1D_eta(sp, 100, 10, -1)
tock()

tick()
for i = 1:1000
    local u, w    
    waveAiry1D_eta(sp, 10, 0.1, -0.1)
    waveAiry1D_eta(sp, 100, 10, -1)
end
tock()

tick()
t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)
tock()

tick()
res = [ waveAiry1D_pPos(sp, ti, 0.1, -0.1 ) for ti in t ]
η = [x[1] for x in res]
px = [x[2] for x in res]
pz = [x[3] for x in res]
# η, px, pz = waveAiry1D_pPos(sp, t, 0.1, -0.1 )
# @show res
tock()

plot(t, η, dpi=330)
savefig(pltFolder*"ts_elevation.png")

plot(t, ϕ, dpi=330)
savefig(pltFolder*"ts_phi.png")

plot(t, u, dpi=330)
savefig(pltFolder*"ts_velu.png")

plot(t, w, dpi=330)
savefig(pltFolder*"ts_velw.png")

plot(t, px, dpi=330)
savefig(pltFolder*"ts_px.png")

plot(t, pz, dpi=330)
savefig(pltFolder*"ts_pz.png")


# end