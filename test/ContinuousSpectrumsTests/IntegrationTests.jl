module IntegrationTests

using Test
using WaveSpec.ContinuousSpectrums

export prepare_test_models

"""
    prepare_test_models()
Prepares a list of implemented continuous spectral models, with tuples: (model_name::String, model_instance::AbstractSpectrum)
"""
function prepare_test_models()
    # Target Parameters
    Hs_target = 3.0
    Tp_target = 2.0
    h_target  = 30.0
    U10_target = 10.0

    Hs_target_2 = 1.0
    Tp_target_2 = 4.0
    λ1 = 2
    λ2 = 3

    # Prepare Models in list of tuples (name, structure)
    models = [("JONSWAP", JONSWAP(Hs_target, Tp_target)), 
              ("TMA", TMA(Hs_target, Tp_target, h_target)),
              ("Donelan", Donelan(Hs_target, Tp_target, U10_target)), 
              ("Bretschneider", Bretschneider(Hs_target, Tp_target)),
              ("OchiHubble", OchiHubble(Hs_target, Tp_target, λ1, Hs_target_2, Tp_target_2, λ2))]

    return models
end

"""
    test_spectral_integration(spec::AbstractSpectrum; rtol=1e-2)

Generic test suite to verify that a spectral model correctly conserves energy
and respects the relationship between Hs and m0.
"""
function test_spectral_integration(spec::AbstractSpectrum; rtol=1e-2)
    # 1. Expected Energy
    Hs_target = get_Hs(spec)
    m0_exact = (Hs_target / 4.0)^2

    # 2. Numerical Integration
    # Ensure your 'integrate' function is dispatched for AbstractSpectrum
    m0_num = integrate(spec, npoints=1000, order=4)
    Hs_calc = 4 * sqrt(m0_num)

    # 3. Validation
    @test Hs_calc ≈ Hs_target rtol = rtol
    @test m0_num ≈ m0_exact rtol = rtol

    # 4. Tail Coverage Test
    # Verify that the default integration range captures the physical energy
    f_min = get_fmin(spec)
    f_max_wide = get_fmax(spec, multiplier=20.0)
    m0_ext = integrate(spec, f_min, f_max_wide, npoints=1000, order=4)
    
    energy_ratio = m0_num / m0_ext
    @test energy_ratio >= 0.99
end


# Prepare models
models = prepare_test_models()
# Test models
for (name, model) in models
    @testset "Spectrum: $name" begin    
        test_spectral_integration(model)
    end
end

end