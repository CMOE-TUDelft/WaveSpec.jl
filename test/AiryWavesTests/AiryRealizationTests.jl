using Test
using WaveSpec

@testset "AiryRealization" begin
    spec = WaveSpec.ContinuousSpectrums.JONSWAP(1.0, 8.0)
    ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
        spec, WaveSpec.SpectralSampling.UniformSampling(),
        0.05, 0.5, 11; mess=false,
    )
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    state = WaveSpec.AiryWaves.AiryState(ds, spread, 50.0)
    realization = WaveSpec.realize(state)

    @test realization isa WaveSpec.AiryRealization
    @test length(realization) == state.nω * state.nθ
    @test realization.h == state.h

    c = realization.components
    @test length(c) == state.nω * state.nθ
    @test c.ω == repeat(state.ω, inner=state.nθ)
    @test c.kx == repeat(state.k, inner=state.nθ) .* repeat(cos.(state.θ), outer=state.nω)
    @test c.ky == repeat(state.k, inner=state.nθ) .* repeat(sin.(state.θ), outer=state.nω)
    @test c.amplitude == vec(WaveSpec.AiryWaves.get_amplitudes(state))
    @test c.phase == vec(WaveSpec.AiryWaves.get_random_phases(state))

    realization2 = WaveSpec.realize(state)
    @test realization2.components.phase == c.phase
    @test realization2.components.kx == c.kx
    @test realization2.components.ky == c.ky

    x, y = 2.3, -0.7
    for i in 1:state.nω, j in 1:state.nθ
        n = (i - 1) * state.nθ + j
        @test c.kx[n] * x + c.ky[n] * y ≈
            state.k[i] * (x * cos(state.θ[j]) + y * sin(state.θ[j]))
    end
end
