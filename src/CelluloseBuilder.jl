module CelluloseBuilder

## exported .jl functions
export cellulosebuilder
export gettingBasisVectors, gettingPBC, unitcell2cartesian, atomsvecString, _XY_trimming_coords, _Z_propagation_coords
export _exporting_XYZfile, _exporting_PDBfile, transformingPBC, _XYZfragments_2_PDB, _cleaning_PDBfragment
export generate_cellulose_topology

# main cellulose builder function
include("./cellulose.jl")

## editing the default atomnames
include("./editing_atomnames.jl")

## crystalographic tools to deal with the crystalline cellulose
include("./crystaltoolkit.jl")

include("./formatter.jl")
include("./exporting_systems.jl")
include("./PBC.jl")

## picking the number of structure fragments inside a VMD structure loading output
include("./picking_fragments.jl")

## cleaning the temporary files
include("./cleaning_tmpfiles.jl")

include("./topology.jl")
include("./vdw-surface.jl")

end
