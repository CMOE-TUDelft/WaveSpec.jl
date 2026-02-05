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
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/CMOE-TUDelft/WaveSpec.jl",
    devbranch="main",
)
