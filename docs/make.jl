using WaveSpec
using Documenter

DocMeta.setdocmeta!(WaveSpec, :DocTestSetup, :(using WaveSpec); recursive=true)

makedocs(;
    modules=[WaveSpec],
    authors="Shagun Agarwal, Pau Manyer, Oriol Colomes",
    repo="https://github.com/CMOE-TUDelft/WaveSpec.jl/blob/{commit}{path}#{line}",
    sitename="WaveSpec.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CMOE-TUDelft.github.io/WaveSpec.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => [
            "Installation" => "guide/installation.md",
            "Getting Started" => "guide/getting_started.md",
             "Tutorials" => "guide/tutorials.md",
        ],
        "API Reference" => [
            "Continuous Spectrums" => "api/continuous_spectrums.md",
            "Spectral Sampling" => "api/spectral_sampling.md",
            "Spectral Spreading" => "api/spectral_spreading.md",
            "Angular Spreading" => "api/angular_spreading.md",
            "Airy Waves" => "api/airy_waves.md",
            "Utils" => "api/utils.md",
        ],
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/CMOE-TUDelft/WaveSpec.jl",
    devbranch="main",
)
