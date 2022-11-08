using CirculantEmbedding
using Documenter

DocMeta.setdocmeta!(CirculantEmbedding, :DocTestSetup, :(using CirculantEmbedding); recursive=true)

makedocs(;
    modules=[CirculantEmbedding],
    authors="Jake Grainger",
    repo="https://github.com/JakeGrainger/CirculantEmbedding.jl/blob/{commit}{path}#{line}",
    sitename="CirculantEmbedding.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JakeGrainger.github.io/CirculantEmbedding.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JakeGrainger/CirculantEmbedding.jl",
    devbranch="main",
)
