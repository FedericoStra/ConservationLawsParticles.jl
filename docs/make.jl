using ConservationLawsParticles
using Documenter

DocMeta.setdocmeta!(ConservationLawsParticles, :DocTestSetup, :(using ConservationLawsParticles); recursive=true)

makedocs(;
    modules=[ConservationLawsParticles],
    authors="Federico Stra <stra.federico@gmail.com> and contributors",
    repo="https://github.com/FedericoStra/ConservationLawsParticles.jl/blob/{commit}{path}#{line}",
    sitename="ConservationLawsParticles.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://FedericoStra.github.io/ConservationLawsParticles.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "examples.md",
    ],
)

deploydocs(;
    repo="github.com/FedericoStra/ConservationLawsParticles.jl",
)
