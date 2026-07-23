using Test
using WaveSpec

@testset "AiryWaves: generate_interpolable_sea" begin
    # Construct a simple AiryState using a RegularWave (monochromatic)
    spec = WaveSpec.ContinuousSpectrums.RegularWave(1.0, 5.0)
    ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec)
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    as = WaveSpec.AiryWaves.AiryState(ds, spread, 50.0)

    # Small evaluation grid
    x = [0.0, 1.0, 2.0]
    y = [0.0, 0.5]
    z = [-1.0, 0.0]
    t = [0.0, 0.5, 1.0]

    # Compute reference grid values with generate_sea
    sea = WaveSpec.AiryWaves.generate_sea(as, x, y, z, t; vars=[:η])
    @test haskey(sea, :η)
    ηgrid = sea[:η]
    # Expand singleton dimensions in the reference grid so it matches the evaluation axes
    target_lens = (length(x), length(y), length(z), length(t))
    cur_shape = size(ηgrid)
    if length(cur_shape) < 4
        newshape = Tuple(vcat(collect(cur_shape), ones(Int, 4 - length(cur_shape))))
        ηgrid = reshape(ηgrid, newshape)
        cur_shape = size(ηgrid)
    end
    reps = ntuple(i -> cur_shape[i] == target_lens[i] ? 1 : target_lens[i], 4)
    if any(r -> r != 1, reps)
        ηgrid = repeat(ηgrid, reps...)
    end

    # Build interpolants with different kernels
    itps_lin = WaveSpec.AiryWaves.generate_interpolable_sea(as, x, y, z, t; vars=[:η], interp=:linear)
    itps_near = WaveSpec.AiryWaves.generate_interpolable_sea(as, x, y, z, t; vars=[:η], interp=:nearest)

    @test haskey(itps_lin, :η) && haskey(itps_near, :η) 

    # Interpolants should reproduce the grid values at the node points
    for ix in eachindex(x), iy in eachindex(y), iz in eachindex(z), it in eachindex(t)
        xv = x[ix]; yv = y[iy]; zv = z[iz]; tv = t[it]
        expected = ηgrid[ix,iy,iz,it]
        @test isapprox(itps_lin[:η](xv,yv,zv,tv), expected; atol=1e-12, rtol=0)
        @test isapprox(itps_near[:η](xv,yv,zv,tv), expected; atol=1e-12, rtol=0)
    end
end

@testset "AiryWaves: change_seed!" begin
    # Regression for the state.spec/state.spectrum field bug: change_seed! must
    # return a working AiryState with the requested seed.
    spec = WaveSpec.ContinuousSpectrums.JONSWAP(1.0, 8.0)
    ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
             spec, WaveSpec.SpectralSampling.UniformSampling(), 0.05, 0.5, 11; mess=false)
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    as = WaveSpec.AiryWaves.AiryState(ds, spread, 50.0)

    as1 = WaveSpec.AiryWaves.change_seed!(as, 12345)
    @test as1 isa WaveSpec.AiryWaves.AiryState
    @test as1.seed == 12345
    # spectrum metadata is carried over unchanged
    @test as1.nω == as.nω && as1.ω == as.ω && as1.k == as.k && as1.h == as.h

    # phases: deterministic per seed, different across seeds
    as2 = WaveSpec.AiryWaves.change_seed!(as, 12345)
    as3 = WaveSpec.AiryWaves.change_seed!(as, 54321)
    @test WaveSpec.AiryWaves.get_random_phases(as1) == WaveSpec.AiryWaves.get_random_phases(as2)
    @test WaveSpec.AiryWaves.get_random_phases(as1) != WaveSpec.AiryWaves.get_random_phases(as3)

    # negative seeds are rejected
    @test_throws ArgumentError WaveSpec.AiryWaves.change_seed!(as, -1)

    # the no-argument variant draws a valid random seed
    as4 = WaveSpec.AiryWaves.change_seed!(as)
    @test as4 isa WaveSpec.AiryWaves.AiryState && as4.seed >= 0
end
