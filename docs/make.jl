using WaveSpec
using Documenter

DocMeta.setdocmeta!(WaveSpec, :DocTestSetup, :(using WaveSpec); recursive=true)

makedocs(;
    modules=[WaveSpec],
    authors="Shagun Agarwal shagun.1994@gmail.com, Thomas Kluwer, Oriol Colomes",
    repo="https://github.com/shagunTUD/WaveSpec.jl/blob/{commit}{path}#{line}",
    sitename="WaveSpec.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://shagunTUD.github.io/WaveSpec.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/shagunTUD/WaveSpec.jl",
    devbranch="main",
)
