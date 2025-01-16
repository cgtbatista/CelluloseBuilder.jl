module CelluloseBuilder

## structures
export XYZ
export XYZs

## functions
export cellulosebuilder
export lattice2basis, PBC, fractional2cartesian, _trimming_xy, _expanding_z
export translate
export writeXYZ, rawXYZ, _exporting_PDBfile, pbcXYZ, xyz2pdb, cleanPDB
export generate_cellulose_topology, get_crystallographic_info

## crystalographic tools to deal with cellulose
include("./crystaltoolkit.jl")

## editing the default atomnames
include("./atomnames.jl")                                       #done
include("./coords.jl")                                          #done

include("./exporting_systems.jl")

## PBC tools
include("./pbc.jl")                                             #done
include("./pbc-misc.jl")

## picking the number of structure fragments inside a VMD structure loading output
include("./fragments.jl")

## cleaning the temporary files
include("./cleaning_tmpfiles.jl")
include("./pdb.jl")
include("./topology.jl")
include("./vdw-surface.jl")

# main cellulose builder function
include("./cellulose.jl")
include("./operators.jl")
include("./utils.jl")

end
