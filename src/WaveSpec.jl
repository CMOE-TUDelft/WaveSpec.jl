module WaveSpec

# Include the foundational utilities first
include("Utils/PhysicalConstants.jl")
include("Utils/Truncation.jl")
include("Utils/Integration.jl")

# Include Continuous Spectrums module
include("ContinuousSpectrums/ContinuousSpectrums.jl")

# Inlcude Spectral Samplings strategies
include("FrequencySampling.jl")

# Include the Discrete Spectrums module
include("SpectralSpreading.jl")

# Include Angular Spreading module
include("AngularSpreading/AngularSpreading.jl")

# Include the Physics/Synthesis module 
include("AiryWaves.jl")


"""
    WaveSpec
A Julia package for generating stochastic sea states using Airy Wave theory,
customizable spectral shapes (JONSWAP, Pierson-Moskowitz), and angular spreading models.
"""
function __init__()
    @info "WaveSpec initialized: Ready for spectral sea state synthesis."
end

end