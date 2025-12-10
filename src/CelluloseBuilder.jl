module CelluloseBuilder

import PDBTools
import LinearAlgebra
import Measurements: Measurement, measurement, ±
import StaticArrays: SVector, @SVector

## structures
export CrystalXYZ, UnitCell

## functions
export cellulose
export lattice2basis, PBC, fractional2cartesian, xy_pruning, z_expansion
export translate
export writeXYZ, rawXYZ, writePDB, pbcXYZ, fragPDBs, cleanPDB
export generate_cellulose_topology, get_crystallographic_info

export surface

export unitcell, crystal, picking_fragments, cleanPDB!, isinverted
export execVMD, execPSFGEN

export polymorph

#####
export atomtypesPSF!, chargesPSF!, atomtypes_generator
export q_cm5, orca_hirshfeld_charges, adjusting_charges

## crystalographic tools to deal with cellulose
include("./crystaltoolkit.jl")

## editing the default atomnames
include("./coords.jl")

## PBC tools
include("./pbc.jl")

## picking the number of structure fragments inside a VMD structure loading output

## cleaning the temporary files
include("./pdb.jl")

# file generators
include("./xyz.jl")

# main cellulose builder function
include("./cellulose.jl")
include("./utils.jl")
include("./topology.jl")        #done without tests
include("./VMD.jl")             #done

include("./opls.jl")

end
