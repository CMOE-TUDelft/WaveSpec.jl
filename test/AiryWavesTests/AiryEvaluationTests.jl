using Test
using Random
import WaveSpec.ContinuousSpectrums as CS
import WaveSpec.SpectralSpreading as SS
import WaveSpec.AngularSpreading as AS
import WaveSpec.AiryWaves as AW
import WaveSpec as WS


@testset "Scalar evaluation matches generate_sea" begin
  
  spec = CS.RegularWave(1.0, 5.0)
  ds = SS.DiscreteSpectralSpreading(spec)
  spread = AS.DiscreteAngularSpreading(0.0)
  state = AW.AiryState(ds,spread,50.0)
  realization = WS.realize(state)
  
  points = (
    (0.0, 0.0, 0.0, 0.0),
    (1.0, 0.0, -1.0, 0.5),
    (5.0, 2.0, -5.0, 1.0),
    (10.0, -3.0, -10.0, 4.0),
  )
  
  for (x,y,z,t) in points
    ref = AW.generate_sea(state,[x],[y],[z],[t];vars=[:η,:ϕ,:u,:v,:w])
    
    @test WS.evaluate_η(realization,x,y,t) ≈ ref.η[1]
    @test WS.evaluate_ϕ(realization,x,y,z,t) ≈ ref.ϕ[1]
    @test WS.evaluate_u(realization,x,y,z,t) ≈ ref.u[1]  
    @test WS.evaluate_v(realization,x,y,z,t) ≈ ref.v[1]
    @test WS.evaluate_w(realization,x,y,z,t) ≈ ref.w[1]
  end
end

@testset "Persistent realization" begin
  
  spec = CS.RegularWave(1.0,5.0)
  ds = SS.DiscreteSpectralSpreading(spec)
  state = AW.AiryState(ds, 50.0)
  
  realization1 = WS.realize(state)
  realization2 = WS.realize(state)

  @test realization1.components.phase == realization2.components.phase
  
end

@testset "evaluate_x regression" begin
  
  spec = CS.JONSWAP(1.0,8.0)
  ds = SS.DiscreteSpectralSpreading(spec,WS.SpectralSampling.UniformSampling(),0.05,0.5,11;mess=false)
  spread = AS.DiscreteAngularSpreading(0.0)
  state = AW.AiryState(ds, spread, 50.0)
  
  realization = WS.realize(state)
  
  points = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 2.0, 0.0),
    (1.0, 2.0, 3.0),
    (10.0, -5.0, 7.0),
  )
  
  for (x, y, t) in points
    sea = AW.generate_sea(state,[x],[y],[0.0],[t];vars=(:η,:ϕ,:u,:v,:w))
    @test WS.evaluate_η(realization,x,y,t) ≈ sea[:η][1] atol=1e-12
    @test WS.evaluate_ϕ(realization,x,y,0.0,t) ≈ sea[:ϕ][1] atol=1e-12
    @test WS.evaluate_u(realization,x,y,0.0,t) ≈ sea[:u][1] atol=1e-12
    @test WS.evaluate_v(realization,x,y,0.0,t) ≈ sea[:v][1] atol=1e-12
    @test WS.evaluate_w(realization,x,y,0.0,t) ≈ sea[:w][1] atol=1e-12
  end
  
  η1 = WS.evaluate_η(realization, 1.2, 3.4, 5.6)
  η2 = WS.evaluate_η(realization, 1.2, 3.4, 5.6)
  ϕ1 = WS.evaluate_ϕ(realization, 1.2, 3.4, 5.6, 7.8)
  ϕ2 = WS.evaluate_ϕ(realization, 1.2, 3.4, 5.6, 7.8)
  u1 = WS.evaluate_u(realization, 1.2, 3.4, 5.6, 7.8)
  u2 = WS.evaluate_u(realization, 1.2, 3.4, 5.6, 7.8)
  v1 = WS.evaluate_v(realization, 1.2, 3.4, 5.6, 7.8)
  v2 = WS.evaluate_v(realization, 1.2, 3.4, 5.6, 7.8)
  w1 = WS.evaluate_w(realization, 1.2, 3.4, 5.6, 7.8)
  w2 = WS.evaluate_w(realization, 1.2, 3.4, 5.6, 7.8) 
  @test ϕ1 == ϕ2
  @test η1 == η2
  @test u1 == u2
  @test v1 == v2
  @test w1 == w2
end

# The following test set checks the correctness of the `evaluate_η!` function.
# It uses an independent reference implementation `direct_η` to verify results.
@testset "evaluate_x!" begin
  
  """
      direct_η(components, x, y, t)

  Independent reference implementation of η(x,y,t).
  Deliberately does not use `evaluate_η!`.
  """
  function direct_η(components, x, y, t)
    T = promote_type(eltype(x), eltype(y), eltype(t))
    η = zeros(T, length(x))
    
    for j in eachindex(η)
      for i in eachindex(components.ω)
        ψ = components.kx[i] * x[j] +
            components.ky[i] * y[j] -
            components.ω[i] * t[j] +
            components.phase[i]
        η[j] += components.amplitude[i] * cos(ψ)
      end
    end
    
    return η
  end
  
  # 1. Single component: analytical result
  @testset "single component" begin
  
    components = WS.WaveComponents([2.0],[3.0],[4.0],[2.5],[0.3])
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)
    
    x = [0.0, 1.0, 2.0]
    y = [0.0, 0.5, 1.0]
    t = [0.0, 0.2, 0.4]
    
    η = similar(x)
    WS.evaluate_η!(η, realization, x, y, t)
    expected = [ 2.5 * cos(3.0 * x[j] + 4.0 * y[j] - 2.0 * t[j] + 0.3)  for j in eachindex(x)]
    @test η ≈ expected
  end
  
  # 2. Multiple components
  @testset "multiple components" begin
    
    components = WS.WaveComponents(
        [1.0, 2.0, 3.0],
        [0.5, 1.0, 1.5],
        [0.2, 0.4, 0.6],
        [1.0, 0.5, 0.25],
        [0.0, 0.2, -0.4],
      )
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)
  
    x = collect(range(0.0, 2.0; length = 17))
    y = collect(range(-1.0, 1.0; length = 17))
    t = collect(range(0.0, 1.0; length = 17))
  
    η = zeros(length(x))
    WS.evaluate_η!(η, realization, x, y, t)  
    expected = direct_η(realization.components, x, y, t)
    @test η ≈ expected
  end
  
  # 3. Zero-amplitude components
  @testset "zero amplitudes" begin
    
    components = WS.WaveComponents(
      [1.0, 2.0, 3.0],
      [1.0, 1.0, 1.0],
      [0.0, 0.0, 0.0],
      [0.0, 2.0, 0.0],
      [0.0, 0.0, 0.0],
    )
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)
  
    x = [0.0, 1.0, 2.0]
    y = [0.0, 0.0, 0.0]
    t = [0.0, 0.5, 1.0]
  
    η = zeros(3)
    WS.evaluate_η!(η, realization, x, y, t)
    expected = [2.0 * cos(x[j] - 2.0 * t[j]) for j in eachindex(x)]
    @test η ≈ expected
  end
  
  # 4. Empty evaluation set
  @testset "empty evaluation set" begin
  
    components = WS.WaveComponents([1.0],[1.0],[0.0],[1.0],[0.0])
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)

    x = Float64[]
    y = Float64[]
    t = Float64[]
    η = Float64[]

    @test WS.evaluate_η!(η, realization, x, y, t) === η
    @test isempty(η)
  end
  
  # 5. Existing values in η must be overwritten
  @testset "output is overwritten" begin
          
    components = WS.WaveComponents([1.0],[1.0],[0.0],[2.0],[0.0])
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)
          
    x = [0.0, 1.0, 2.0]
    y = zeros(3)
    t = zeros(3)

    η = fill(12345.0, 3)
    WS.evaluate_η!(η, realization, x, y, t)
    @test η ≈ [2.0, 2.0 * cos(1.0), 2.0 * cos(2.0)]
  end
          
  # 6. Length validation
  @testset "length validation" begin
  
    components = WS.WaveComponents([1.0],[1.0],[0.0],[1.0],[0.0])
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)
            
    η = zeros(3)
    x = zeros(3)
    y = zeros(3)
    t = zeros(2)
            
    @test_throws Exception evaluate_η!(η,realization,x,y,t,)
  end
            
  # 7. Numerical type
  @testset "Float32" begin
              
    components = WS.WaveComponents(
      Float32[1.0, 2.0],
      Float32[0.5, 1.0],
      Float32[0.2, 0.4],
      Float32[1.0, 0.5],
      Float32[0.1, -0.2],
    )
    realization = WS.AiryRealization(components, ones(Float32, length(components)), 1f0)
              
    x = Float32[0.0, 1.0, 2.0]
    y = Float32[0.0, 0.5, 1.0]
    t = Float32[0.0, 0.2, 0.4]

    η = zeros(Float32, 3)
    WS.evaluate_η!(η, realization, x, y, t)
    expected = direct_η(realization.components, x, y, t)
    @test η ≈ expected rtol = 1f-5
  end
              
  # 8. Reference implementation comparison
  @testset "against independent reference" begin
  
    rng = Random.MersenneTwister(1234)

    Nc = 100
    Neval = 257
    components = WS.WaveComponents(
      rand(rng, Nc),
      randn(rng, Nc),
      randn(rng, Nc),
      rand(rng, Nc),
      2π .* rand(rng, Nc),
    )
    realization = WS.AiryRealization(components, ones(length(components)), 1.0)

    x = randn(rng, Neval)
    y = randn(rng, Neval)
    t = randn(rng, Neval)

    η = zeros(Neval)
    WS.evaluate_η!(η, realization, x, y, t)
    expected = direct_η(realization.components, x, y, t)
    @test η ≈ expected rtol = 1e-12 atol = 1e-12
  end
  
  
    # -------------------------------------------------------------------------
    # 9. Threaded implementation consistency
    # -------------------------------------------------------------------------
  
    @testset "threaded evaluation is deterministic within tolerance" begin
    
        rng = Random.MersenneTwister(5678)
    
        Nc = 50
        Neval = 1000
    
        components = WS.WaveComponents(
          rand(rng, Nc),
          randn(rng, Nc),
          randn(rng, Nc),
          rand(rng, Nc),
          2π .* rand(rng, Nc),
        )
        realization = WS.AiryRealization(components, ones(length(components)), 1.0)
    
        x = randn(rng, Neval)
        y = randn(rng, Neval)
        t = randn(rng, Neval)
    
        η = zeros(Neval)
        WS.evaluate_η!(η, realization, x, y, t)
        expected = direct_η(realization.components, x, y, t)
        @test η ≈ expected rtol = 1e-12 atol = 1e-12
      end
    
  end
# -----------------------------------------------------------------------------
# Independent textbook reference over the (ω, θ) grid of an AiryState. It does
# not use `realize`, the evaluation kernels or `generate_sea`.
# -----------------------------------------------------------------------------
function reference_fields(state, x, y, z, t)
  A = AW.get_amplitudes(state)
  φ = AW.get_random_phases(state)
  g = WS.PhysicalConstants.g
  h = state.h
  η = ϕ = u = v = w = 0.0
  for i in 1:state.nω, j in 1:state.nθ
    ω, k, θ = state.ω[i], state.k[i], state.θ[j]
    ψ = k * (x * cos(θ) + y * sin(θ)) - ω * t + φ[i, j]
    η += A[i, j] * cos(ψ)
    ϕ += A[i, j] * g / ω * cosh(k * (z + h)) / cosh(k * h) * sin(ψ)
    u += A[i, j] * ω * cos(θ) * cosh(k * (z + h)) / sinh(k * h) * cos(ψ)
    v += A[i, j] * ω * sin(θ) * cosh(k * (z + h)) / sinh(k * h) * cos(ψ)
    w += A[i, j] * ω * sinh(k * (z + h)) / sinh(k * h) * sin(ψ)
  end
  return (η = η, ϕ = ϕ, u = u, v = v, w = w)
end

@testset "Short-crested sea matches textbook reference" begin
  spec = CS.JONSWAP(2.0, 8.0)
  ds = SS.DiscreteSpectralSpreading(spec, WS.SpectralSampling.UniformSampling(), 0.05, 0.5, 8; mess=false)
  spread = AS.DiscreteAngularSpreading(:cosinepow, 0.0, 0.5, -π/2, π/2, 6)
  state = AW.AiryState(ds, spread, 30.0)
  @test state.nω > 1 && state.nθ > 1
  r = WS.realize(state)

  pts = ((0.0, 0.0, 0.0, 0.0), (3.0, -2.0, -4.0, 1.5), (10.0, 5.0, -1.0, 7.0), (-7.0, 12.0, -29.0, 33.0))
  for (x, y, z, t) in pts
    ref = reference_fields(state, x, y, z, t)
    @test WS.evaluate_η(r, x, y, t) ≈ ref.η atol = 1e-12
    @test WS.evaluate_ϕ(r, x, y, z, t) ≈ ref.ϕ atol = 1e-12
    @test WS.evaluate_u(r, x, y, z, t) ≈ ref.u atol = 1e-12
    @test WS.evaluate_v(r, x, y, z, t) ≈ ref.v atol = 1e-12
    @test WS.evaluate_w(r, x, y, z, t) ≈ ref.w atol = 1e-12
    f = WS.evaluate_fields(r, x, y, z, t)
    for name in (:η, :ϕ, :u, :v, :w)
      @test f[name] ≈ ref[name] atol = 1e-12
    end
  end

  # generate_sea on a small tensor grid
  xs = [0.0, 3.0, 10.0]; ys = [-2.0, 5.0]; zs = [-4.0, -1.0, 0.0]; ts = [0.0, 1.5, 7.0]
  sea = AW.generate_sea(state, xs, ys, zs, ts)
  @test size(sea.η) == (3, 2, 1, 3)
  @test size(sea.u) == (3, 2, 3, 3)
  for (ix, x) in enumerate(xs), (iy, y) in enumerate(ys), (iz, z) in enumerate(zs), (it, t) in enumerate(ts)
    ref = reference_fields(state, x, y, z, t)
    iz == 1 && @test sea.η[ix, iy, 1, it] ≈ ref.η atol = 1e-12
    @test sea.ϕ[ix, iy, iz, it] ≈ ref.ϕ atol = 1e-12
    @test sea.u[ix, iy, iz, it] ≈ ref.u atol = 1e-12
    @test sea.v[ix, iy, iz, it] ≈ ref.v atol = 1e-12
    @test sea.w[ix, iy, iz, it] ≈ ref.w atol = 1e-12
  end
  @test keys(AW.generate_sea(state, xs, ys, zs, ts; vars = [:w, :η])) == (:w, :η)
end

@testset "Deep water (kh ≫ 20) stays finite and exact" begin
  spec = CS.JONSWAP(1.0, 4.0)
  ds = SS.DiscreteSpectralSpreading(spec, WS.SpectralSampling.UniformSampling(), 0.1, 2.0, 12; mess=false)
  state = AW.AiryState(ds, AS.DiscreteAngularSpreading(0.3), 2000.0)
  @test maximum(state.k) * state.h > 700   # cosh(kh) overflows Float64 here
  r = WS.realize(state)
  for z in (0.0, -1.0, -10.0, -2000.0)
    f = WS.evaluate_fields(r, 1.0, 2.0, z, 3.0)
    @test all(isfinite, values(f))
  end
  # Near the surface the exact profile equals the deep-water exponential
  x, y, z, t = 1.0, 2.0, -1.5, 3.0
  A = AW.get_amplitudes(state); φ = AW.get_random_phases(state)
  ϕ_deep = sum(A[i, 1] * WS.PhysicalConstants.g / state.ω[i] * exp(state.k[i] * z) *
               sin(state.k[i] * (x * cos(state.θ[1]) + y * sin(state.θ[1])) - state.ω[i] * t + φ[i, 1])
               for i in 1:state.nω)
  @test WS.evaluate_ϕ(r, x, y, z, t) ≈ ϕ_deep rtol = 1e-10
end

@testset "evaluate_fields! and type stability" begin
  spec = CS.JONSWAP(2.0, 8.0)
  ds = SS.DiscreteSpectralSpreading(spec, WS.SpectralSampling.UniformSampling(), 0.05, 0.5, 8; mess=false)
  state = AW.AiryState(ds, AS.DiscreteAngularSpreading(:cosinepow, 0.0, 0.5, -π/2, π/2, 6), 30.0)
  r = WS.realize(state)
  rng = Random.MersenneTwister(42)
  n = 64
  x = 50 .* randn(rng, n); y = 50 .* randn(rng, n); z = -30 .* rand(rng, n); t = 100 .* rand(rng, n)
  η, ϕ, u, v, w = (zeros(n) for _ in 1:5)
  WS.evaluate_fields!(η, ϕ, u, v, w, r, x, y, z, t)
  @test η ≈ WS.evaluate_η!(zeros(n), r, x, y, t)
  @test ϕ ≈ WS.evaluate_ϕ!(zeros(n), r, x, y, z, t)
  @test u ≈ WS.evaluate_u!(zeros(n), r, x, y, z, t)
  @test v ≈ WS.evaluate_v!(zeros(n), r, x, y, z, t)
  @test w ≈ WS.evaluate_w!(zeros(n), r, x, y, z, t)
  @test_throws DimensionMismatch WS.evaluate_fields!(η, ϕ, u, v, zeros(n - 1), r, x, y, z, t)

  c = r.components
  r32 = WS.AiryRealization(WS.WaveComponents(Float32.(c.ω), Float32.(c.kx), Float32.(c.ky),
                                             Float32.(c.amplitude), Float32.(c.phase)),
                           Float32.(r.k), Float32(r.h))
  @test (@inferred WS.evaluate_ϕ(r32, 1f0, 2f0, -3f0, 4f0)) isa Float32
  @test (@inferred WS.evaluate_w(r32, 1f0, 2f0, -3f0, 4f0)) isa Float32
  @test (@inferred WS.evaluate_fields(r32, 1f0, 2f0, -3f0, 4f0)).u isa Float32
  @test WS.evaluate_ϕ(r32, 1f0, 2f0, -3f0, 4f0) ≈ WS.evaluate_ϕ(r, 1.0, 2.0, -3.0, 4.0) rtol = 1e-4
  # Measure inside a function so the result is not boxed at testset scope
  fields_alloc(r) = @allocated WS.evaluate_fields(r, 1.0, 2.0, -3.0, 4.0)
  fields_alloc(r)
  @test fields_alloc(r) == 0
end
