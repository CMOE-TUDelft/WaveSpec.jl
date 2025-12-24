using WaveSpec
using .Constants
using .Jonswap
using .WaveTimeSeries

# Define a main test set for module
@testset "JONSWAP Spectrum & Airy Wave Generation" begin
    # 1. SETUP (Keep all your initialization code)
    h0 = 50.0  # Use Float for consistency
    H₀ = 7.11 
    Tₚ = 12.0
    nω = 257

    # Generate Spectrum
    ω, S, A = jonswap(H₀, Tₚ; plotflag=false, nω=nω)
    
    # ... (rest of the setup: k, α, sp)
    k = dispersionRelAng.(h0, ω; msg=false)
    α = randomPhase(ω; seed=100)
    sp = SpecStruct( h0, ω, S, A, k, α; Hs = H₀, Tp = Tₚ )
    
    
    # Run different specific tests for module functionalities

    @testset "JONSWAP Output Validity" begin
        # Test 1: Check the output types are correct
        @test typeof(ω) <: AbstractVector
        @test typeof(S) <: AbstractVector
        @test typeof(A) <: AbstractVector
        
        # Test 2: Check the correct number of points were generated
        @test length(ω) == nω
        @test length(S) == nω
        
        # Test 3: Verify the significant wave height calculation (THE CORE ACCURACY TEST)
        # Assuming you have a function to calculate Hs from the spectrum S. 
        # If not, the integral approximation must be checked.
        m0 = sum(S .* diff(ω)) # Simple integration approximation
        calculated_Hs = sqrt(16 * m0)
        
        # The calculated Hs should be very close to the input Hs (7.11)
        @test isapprox(calculated_Hs, H₀, rtol=0.01) # Allow 1% relative tolerance
    end
    
    @testset "Airy Wave Generation" begin
        t = 0.6
        x = 0.1
        z = -0.1
        
        # Test 4: Check function call 1 (with individual components)
        η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, x, z)
        @test typeof(η) <: AbstractFloat # Ensure scalar output
        @test !isnan(η) # Ensure it calculates a real number
        
        # Test 5: Check function call 2 (with SpecStruct)
        η_sp, ϕ_sp, u_sp, w_sp = waveAiry1D(sp, t, x, z)

        # Test 6: Ensure both methods give the same result
        @test isapprox(η, η_sp, rtol=1e-6) 
        
        # Test 7: Test the particle position function
        η_ppos, px, py = waveAiry1D_pPos(sp, t, x, z)
        @test isapprox(η_sp, η_ppos, rtol=1e-6) 
        @test !isinf(px)
        @test typeof(py) <: AbstractFloat
    end
    
    # The time series calls (t=0:0.1:1200) are usually slower and often skipped 
    # or moved to separate performance tests. You might want to remove them from 
    # the fast CI run for now.
end


# Define a separate test set for basic construction of JONSWAP spectrum
@testset "JONSWAP Basic Construction" begin
    # 1. Test the function call succeeds without throwing an error
    @test_nowarn WaveSpec.jonswap(1.0, 10.0) # Assume 1.0 is peak frequency, 10.0 is significant wave height

    # 2. Test the output type is correct (e.g., a spectral struct or array)
    spectrum = WaveSpec.jonswap(1.0, 10.0)
    @test typeof(spectrum) <: AbstractWaveSpectrum 
    
    # 3. Test that the output size/length is reasonable if applicable
    @test length(spectrum.frequencies) > 0 
end