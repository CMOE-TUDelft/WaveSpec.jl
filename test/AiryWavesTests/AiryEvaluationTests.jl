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
    realization = WS.AiryRealization(components, [0.0], 0.0)
    
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
    realization = WS.AiryRealization(components,[0.0],0.0)
  
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
    realization = WS.AiryRealization(components,[0.0],0.0)
  
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
    realization = WS.AiryRealization(components, [0.0], 0.0)

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
    realization = WS.AiryRealization(components, [0.0], 0.0)
          
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
    realization = WS.AiryRealization(components, [0.0], 0.0)
            
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
    realization = WS.AiryRealization(components,Float32[0.0],Float32(0.0))
              
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
    realization = WS.AiryRealization(components,[0.0],0.0)

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
        realization = WS.AiryRealization(components,[0.0],0.0,)
    
        x = randn(rng, Neval)
        y = randn(rng, Neval)
        t = randn(rng, Neval)
    
        η = zeros(Neval)
        WS.evaluate_η!(η, realization, x, y, t)
        expected = direct_η(realization.components, x, y, t)
        @test η ≈ expected rtol = 1e-12 atol = 1e-12
      end
    
  end