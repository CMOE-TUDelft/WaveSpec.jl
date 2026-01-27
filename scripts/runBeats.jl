# Simulate single and two freq to observe the beats effect

module run_wave

using   WaveSpec
using   Plots
using   .Constants
using   .WaveTimeSeries

pltFolder = "./output/"

ω = [2*pi/9.9, 2*pi/10.1]
A = [1,1]
α = [0.0, 0.0]
k = dispersionRelAng.(30, ω)

t=0:0.1:500
#η, ϕ, u, w = waveAiry1D(30, ω, A, k, α, t)
η, ϕ, u, w = waveStokes2nd1D(30, ω[1], A[1], k[1], α[1], t)

@show k

scatter(ω, A)
savefig(pltFolder*"spec.png")

plot(t, η)
savefig(pltFolder*"ts_elevation.png")

plot(t, ϕ)
savefig(pltFolder*"ts_phi.png")

plot(t, u)
savefig(pltFolder*"ts_velu.png")

plot(t, w)
savefig(pltFolder*"ts_velw.png")

end