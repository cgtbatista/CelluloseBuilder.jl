module CelluloseBuilder

## functions
export cellulosebuilder
export gettingBasisVectors, gettingPBC, fractional2cartesian, atomnames, _trimming_xy, _expanding_z
export transformASU
export writeXYZ, rawXYZ, _exporting_PDBfile, transformingPBC, xyz2pdb, cleanPDB
export generate_cellulose_topology, get_crystallographic_info

## structures
struct XYZ
    atoms::Vector{String}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct XYZs
    atoms::Vector{Vector{String}}
    x::Vector{Vector{Float64}}
    y::Vector{Vector{Float64}}
    z::Vector{Vector{Float64}}
end

export XYZ, XYZs

## editing the default atomnames
include("./atomnames.jl")
include("./expanding.jl")
include("./coords.jl")
## crystalographic tools to deal with the crystalline cellulose
include("./crystaltoolkit.jl")

include("./exporting_systems.jl")
include("./PBC.jl")

## picking the number of structure fragments inside a VMD structure loading output
include("./fragments.jl")

## cleaning the temporary files
include("./cleaning_tmpfiles.jl")
include("./pdb.jl")
include("./topology.jl")
include("./vdw-surface.jl")

# main cellulose builder function
include("./cellulose.jl")

include("./utils.jl")

end
