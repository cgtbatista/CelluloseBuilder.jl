module CelluloseBuilder

## structures
export CrystalXYZ, UnitCell

## functions
export cellulosebuilder
export lattice2basis, PBC, fractional2cartesian, xy_pruning, z_expansion
export translate
export writeXYZ, rawXYZ, _exporting_PDBfile, pbcXYZ, fragPDBs, cleanPDB
export generate_cellulose_topology, get_crystallographic_info

export unitcell, crystal, picking_fragments

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

## cleaning the temporary files
include("./cleaning_tmpfiles.jl")
include("./pdb.jl")
include("./topology.jl")
include("./vdw-surface.jl")

# file generators
include("./xyz.jl")

# main cellulose builder function
include("./cellulose.jl")
include("./operators.jl")
include("./utils.jl")

end
