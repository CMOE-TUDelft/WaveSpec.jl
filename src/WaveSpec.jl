module WaveSpec

# Include the foundational utilities first
include("Utils/PhysicalConstants.jl")
include("Utils/Integration.jl")
include("Utils/Truncation.jl")

# Include Continuous Spectrums module
include("ContinuousSpectrums/ContinuousSpectrums.jl")

# Inlcude Spectral Samplings strategies
include("SpectralSampling.jl")

# Include the Discrete Spectrums module
include("SpectralSpreading.jl")

# Include Angular Spreading module
include("AngularSpreading/AngularSpreading.jl")

# Include the wave-component representation
include("WaveComponents.jl")

# Include Airy Waves module
include("AiryWaves.jl")
include("AiryRealization.jl")
include("AiryEvaluation.jl")

# Include signal generation and treatment tools module
include("Utils/Signal.jl")


"""
    WaveSpec
A Julia package for generating stochastic sea states using Airy Wave theory,
customizable spectral shapes (JONSWAP, Pierson-Moskowitz), and angular spreading models.
"""
function __init__()
    @info "WaveSpec initialized: Ready for spectral sea state synthesis."
end

end