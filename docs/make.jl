using CelluloseBuilder
using Documenter

DocMeta.setdocmeta!(CelluloseBuilder, :DocTestSetup, :(using CelluloseBuilder); recursive=true)

makedocs(;
    modules=[CelluloseBuilder],
    authors="cgtbatista <c203748@dac.unicamp.br> and contributors",
    sitename="CelluloseBuilder.jl",
    format=Documenter.HTML(;
        canonical="https://cgtbatista.github.io/CelluloseBuilder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cgtbatista/CelluloseBuilder.jl",
    devbranch="main",
)
