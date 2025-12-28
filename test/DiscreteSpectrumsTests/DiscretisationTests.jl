module DiscretisationTests

using Test
using WaveSpec.ContinuousSpectrums
using WaveSpec.SpectralSpreading


@testset "Discrete Spectrum Bins" begin
    spec_base = JONSWAP(3.0, 10.0)
    # Create a discrete container with 100 frequencies
    # Using the bin-edge logic we implemented (Midpoint/Central)
    ds = DiscreteSpectrum(spec_base, f_min=0.02, f_max=0.5, nf=101)
    
    # Check that we have nf-1 bins
    @test length(get_bin_centers(ds)) == 100
    
    # Ensure bandwidths sum to the total range
    @test sum(get_bandwidths(ds)) ≈ (0.5 - 0.02)
end

end # module