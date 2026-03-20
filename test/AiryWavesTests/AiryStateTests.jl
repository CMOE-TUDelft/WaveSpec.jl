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
    itps_cub = WaveSpec.AiryWaves.generate_interpolable_sea(as, x, y, z, t; vars=[:η], interp=:cubic)

    @test haskey(itps_lin, :η) && haskey(itps_near, :η) && haskey(itps_cub, :η)

    # Interpolants should reproduce the grid values at the node points
    for ix in eachindex(x), iy in eachindex(y), iz in eachindex(z), it in eachindex(t)
        xv = x[ix]; yv = y[iy]; zv = z[iz]; tv = t[it]
        expected = ηgrid[ix,iy,iz,it]
        @test isapprox(itps_lin[:η](xv,yv,zv,tv), expected; atol=1e-12, rtol=0)
        @test isapprox(itps_near[:η](xv,yv,zv,tv), expected; atol=1e-12, rtol=0)
        @test isapprox(itps_cub[:η](xv,yv,zv,tv), expected; atol=1e-8)
    end
end
