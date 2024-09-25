module CelluloseBuilder

import StaticArrays

import PDBTools
import MolSimToolkit
import Packmol

import LinearAlgebra: cross, norm, dot

# using PDBTools, MolSimToolkit, StaticArrays, Packmol, LinearAlgebra, Test
## exported .jl functions

export cellulosebuilder

export gettingBasisVectors, gettingPBC, unitcell2cartesian, atomsvecString, _XY_trimming_coords, _Z_propagation_coords
export _exporting_XYZfile, _exporting_PDBfile, transformingPBC, _XYZfragments_2_PDB, _cleaning_PDBfragment
export generate_cellulose_topology
export updating_segid

# Write your package code here.
##const DEFAULT_CARB_TOPOLOGY_FILE = "$(@__DIR__)/toppar/cellulose.rtf"

# future input format
struct CelluloseBuilderInput
    phase::String
    pbc
    covalent::Bool
    topologyfile::String
    atomnames::Vector{String}
    asym_coords::Vector{Vector{Float64}}
    uc_parameters::Vector{Vector{Float64}}
end
export CelluloseBuilderInput

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

## Misc. tools to deal help the system building
include("./tools-assembly.jl")
include("./tools-maths.jl")
include("./tools-patching.jl")
include("./tools-pdb.jl")
include("./tools-topology.jl")

## getting the surface
include("./vdw-surface.jl")
include("./residue-library.jl")

end