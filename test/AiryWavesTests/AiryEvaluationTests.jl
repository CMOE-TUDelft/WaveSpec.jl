using Test
using WaveSpec

@testset "Scalar evaluation matches generate_sea" begin

    spec = WaveSpec.ContinuousSpectrums.RegularWave(1.0, 5.0)

    ds =
        WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
            spec
        )

    spread =
        WaveSpec.AngularSpreading.DiscreteAngularSpreading(
            0.0
        )

    state =
        WaveSpec.AiryWaves.AiryState(
            ds,
            spread,
            50.0
        )

    realization = WaveSpec.realize(state)

    points = (
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, -1.0, 0.5),
        (5.0, 2.0, -5.0, 1.0),
        (10.0, -3.0, -10.0, 4.0),
    )

    for (x,y,z,t) in points

        ref =
            WaveSpec.AiryWaves.generate_sea(
                state,
                [x],
                [y],
                [z],
                [t];
                vars=[:η,:ϕ,:u,:v,:w]
            )

        @test WaveSpec.evaluate_η(
            realization,
            x,y,t
        ) ≈ ref.η[1]

        @test WaveSpec.evaluate_ϕ(
            realization,
            x,y,z,t
        ) ≈ ref.ϕ[1]

        @test WaveSpec.evaluate_u(
            realization,
            x,y,z,t
        ) ≈ ref.u[1]

        @test WaveSpec.evaluate_v(
            realization,
            x,y,z,t
        ) ≈ ref.v[1]

        @test WaveSpec.evaluate_w(
            realization,
            x,y,z,t
        ) ≈ ref.w[1]

    end
end

@testset "Persistent realization" begin

    spec = WaveSpec.ContinuousSpectrums.RegularWave(
        1.0,
        5.0
    )

    ds =
        WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
            spec
        )

    state =
        WaveSpec.AiryWaves.AiryState(
            ds,
            50.0
        )

    r1 = WaveSpec.realize(state)
    r2 = WaveSpec.realize(state)

    @test r1.components.phase == r2.components.phase

end

@testset "evaluate_η regression" begin

    spec = WaveSpec.ContinuousSpectrums.JONSWAP(
        1.0,
        8.0,
    )

    ds =
        WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
            spec,
            WaveSpec.SpectralSampling.UniformSampling(),
            0.05,
            0.5,
            11;
            mess=false,
        )

    spread =
        WaveSpec.AngularSpreading.DiscreteAngularSpreading(
            0.0,
        )

    state =
        WaveSpec.AiryWaves.AiryState(
            ds,
            spread,
            50.0,
        )

    realization = WaveSpec.realize(state)

    points = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (1.0, 2.0, 3.0),
        (10.0, -5.0, 7.0),
    )

    for (x, y, t) in points

        sea =
            WaveSpec.generate_sea(
                state,
                [x],
                [y],
                [0.0],
                [t];
                vars=(:η,),
            )

        η_old = sea[:η][1]

        η_new =
            WaveSpec.evaluate_η(
                realization,
                x,
                y,
                t,
            )

        @test η_new ≈ η_old atol=1e-12
    end

    η1 =
        WaveSpec.evaluate_η(
            realization,
            1.2,
            3.4,
            5.6,
        )

    η2 =
        WaveSpec.evaluate_η(
            realization,
            1.2,
            3.4,
            5.6,
        )

    @test η1 == η2
end
