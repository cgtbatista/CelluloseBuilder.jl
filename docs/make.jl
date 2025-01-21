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


## To help documentation
## https://ocw.mit.edu/courses/12-108-structure-of-earth-materials-fall-2004/9df654250315660f294bf6c9acd49ae1_lec7.pdf
## https://physics.stackexchange.com/questions/63511/what-is-the-difference-between-lattice-vectors-and-basis-vectors
## http://pd.chem.ucl.ac.uk/pdnn/symm1/trans1.htm
## https://eng.libretexts.org/Bookshelves/Materials_Science/TLP_Library_I/15%3A_Crystallography/15.04%3A_Section_4-
## https://zhanggroup.org/SSIPe/pdb_atom_format.html
