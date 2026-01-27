module run_wave
include("../src/jonswapThomas.jl")

using   WaveSpec
using   Plots
using   .Wave
using   .Constants
using   .Jonswap
using   .WaveTimeSeries

pltFolder = "./output/"

h0 = 50 #still-water depth
H₀ = 7.11 #m
Tₚ = 12 #s

# Thomas Code
nf = 257
fpeak = 1/Tₚ
f = range(0.01, 33*fpeak/2/π, nf)
α, μ, k, kₛ, fₛ, ηₛ, jon = waveparameters(f, fpeak, h0)

Jonswapplot = bar(f, jon)
title!("Jonswap variance density spectrum, fp ="*"$fpeak")
xlabel!("frequency [Hz]")
ylabel!("eta² [m²/Hz]") 
savefig(Jonswapplot, pltFolder*"Jonswap_spectrum_"*"$nf"*"modes_"*"$fpeak"*"fpeak.png")

Amplplot = scatter(f, μ)
title!("Jonswap amplitude spectrum, fp ="*"$fpeak")
xlabel!("frequency [Hz]")
ylabel!("eta [m]") 
savefig(Amplplot, pltFolder*"Jon_amplitude_spectrum_"*"$nf"*"modes_"*"$fpeak"*"fpeak.png")


# Revised
ω, S, A = jonswap(H₀, Tₚ; 
    plotflag=true, plotloc=pltFolder, nω = nf)

k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω; seed=100)

scatter(ω, k)
savefig(pltFolder*"dispersion.png")

scatter(ω, α)
savefig(pltFolder*"phase.png")

t=0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t)

t=0:0.1:0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t)

plot(t, η)
savefig(pltFolder*"ts_elevation.png")

plot(t, ϕ)
savefig(pltFolder*"ts_phi.png")

plot(t, u)
savefig(pltFolder*"ts_velu.png")

plot(t, w)
savefig(pltFolder*"ts_velw.png")

end