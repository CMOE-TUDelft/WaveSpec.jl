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
    @test c.amplitude == vec(permutedims(WaveSpec.AiryWaves.get_amplitudes(state)))
    @test c.phase == vec(permutedims(WaveSpec.AiryWaves.get_random_phases(state)))

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

@testset "AiryRealization component pairing (short-crested)" begin
    spec = WaveSpec.ContinuousSpectrums.JONSWAP(2.0, 8.0)
    ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
        spec, WaveSpec.SpectralSampling.UniformSampling(),
        0.05, 0.5, 6; mess=false,
    )
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(:cosinepow, 0.0, 0.5, -π/2, π/2, 5)
    state = WaveSpec.AiryWaves.AiryState(ds, spread, 30.0)
    @test state.nω > 1 && state.nθ > 1

    c = WaveSpec.realize(state).components
    A = WaveSpec.AiryWaves.get_amplitudes(state)
    φ = WaveSpec.AiryWaves.get_random_phases(state)
    for i in 1:state.nω, j in 1:state.nθ
        n = (i - 1) * state.nθ + j
        @test c.ω[n] == state.ω[i]
        @test c.kx[n] ≈ state.k[i] * cos(state.θ[j])
        @test c.amplitude[n] == A[i, j]
        @test c.phase[n] == φ[i, j]
    end
end

@testset "AiryRealization validates k length" begin
    comps = WaveSpec.WaveComponents([1.0, 2.0], [0.1, 0.4], [0.0, 0.0], [1.0, 1.0], [0.0, 0.0])
    @test_throws DimensionMismatch WaveSpec.AiryRealization(comps, [0.1], 10.0)
end
