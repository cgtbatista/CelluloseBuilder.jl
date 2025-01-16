var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CelluloseBuilder","category":"page"},{"location":"#CelluloseBuilder","page":"Home","title":"CelluloseBuilder","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for CelluloseBuilder.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CelluloseBuilder]","category":"page"},{"location":"#CelluloseBuilder.PBC-Tuple{Int64, Int64, Int64}","page":"Home","title":"CelluloseBuilder.PBC","text":"PBC(xsize::Int64, ysize::Int64, zsize::Int64; phase=\"Iβ\", pbc=nothing)\n\nReturns (ncells, lattice) where ncells is the number of unit cells along x, y and z axes and lattice is the unit cell vector of this cell. This function is used to define the periodic boundary conditions required to build the cellulose crystal.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.PBC-Tuple{Int64, Int64}","page":"Home","title":"CelluloseBuilder.PBC","text":"PBC(monomers::Int64, chains::Int64; phase=\"Iβ\", layer=\"center\")\n\nDefine the crystallographic space to craft a cellulose layer with a degree of β-glucose units (monomers) and a number of cellulose chains (chains).\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.PBC-Tuple{Int64}","page":"Home","title":"CelluloseBuilder.PBC","text":"PBC(monomers::Int64; phase=\"Iβ\")\n\nDefine the crystallographic space to build the cellulose fibril crystal with a degree of β-glucose units (monomers).\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.PBC-Tuple{Vector{Int64}}","page":"Home","title":"CelluloseBuilder.PBC","text":"PBC(xyzsizes::Vector{Int64}; phase=\"Iβ\", pbc=nothing)\n\nSee PBC(xsize::Int64, ysize::Int64, zsize::Int64; phase=\"Iβ\", pbc=nothing).\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.cellulosebuilder-Tuple{Int64}","page":"Home","title":"CelluloseBuilder.cellulosebuilder","text":"cellulosebuilder(monomers::Int64; phase=\"Iβ\", fibril=nothing, ...)\n\nBuilds a cellulose fibril with a given number of cellobiose units based on the number of monomers.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.cleanPDB-Tuple{String, Vector{String}}","page":"Home","title":"CelluloseBuilder.cleanPDB","text":"cleanPDB(pdbname::String, atoms::Vector{String}, new_pdbname::String)\n\nClean the raw PDB file and export a new one with the right configuration needed for CHARMM force field.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.fractional2cartesian-Tuple{Vector{Int64}, String}","page":"Home","title":"CelluloseBuilder.fractional2cartesian","text":"fractional2cartesian(unitcell::Vector{Int64}, phase::String)\n\nConvert the fractional unit cell coordinates to the cartesian coordinates. The return is the x, y, and z cartesian coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.fragPDBs-Tuple{String}","page":"Home","title":"CelluloseBuilder.fragPDBs","text":"fragPDBs(xyzfile::String, selection::String; vmd=\"vmd\")\n\nConvert each fragment on XYZ file on an unique raw PDB file. The filename will be exported, so there is a way to edit the PDB data.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.lattice2basis-Tuple{Vector{Int64}, String}","page":"Home","title":"CelluloseBuilder.lattice2basis","text":"lattice2basis(lattice::Vector{Int64}, phase::String)\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.lattice2basis-Tuple{Vector{Int64}, Vector{Vector{Float64}}}","page":"Home","title":"CelluloseBuilder.lattice2basis","text":"lattice2basis(lattice::Vector{Int64}, parameters::Vector{Vector{Float64}})\n\nConvert the lattice vector to the basis vectors of the unit cell on the cartesian space.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.monolayer-Tuple{Vector{Int64}}","page":"Home","title":"CelluloseBuilder.monolayer","text":"monolayer(xyzsizes::Vector{Int64}; phase=\"Iα\", structure=\"monolayer\")\n\nReturns the fragments selection needed to the monolayer structure.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.pbcXYZ-Tuple{Int64, Int64, Int64}","page":"Home","title":"CelluloseBuilder.pbcXYZ","text":"pbcXYZ(nfrag::Int64, xsize::Int64, ysize::Int64; phase=\"Iβ\", pbc=nothing, xyzfile=\"file.xyz\", vmd=\"vmd\")\n\nApplies the periodic boundary conditions on the cellulose crystal over xy-plane using the number of fragments (e.g. chains) as parameter. It returns cellulose XYZ file needed to prepare the PDB and PSF files.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.pbcdata-Tuple{Int64, Int64, Int64}","page":"Home","title":"CelluloseBuilder.pbcdata","text":"pbcdata(xsize::Int64, ysize::Int64, zsize::Int64)\n\nDictionary useful for (a, b, c) definition of the crystal.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.pbcdata-Tuple{Int64}","page":"Home","title":"CelluloseBuilder.pbcdata","text":"pbcdata(monomers::Int64; chains=nothing)\n\nDictionary useful for cellulose definition of the layers and fibril.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.translate-Tuple{String}","page":"Home","title":"CelluloseBuilder.translate","text":"translate(phase::String)\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.translate-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, String}","page":"Home","title":"CelluloseBuilder.translate","text":"translate(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, phase::String)\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.translate-Tuple{Vector{Vector{Float64}}, String}","page":"Home","title":"CelluloseBuilder.translate","text":"translate(coords::Vector{Vector{Float64}}, phase::String)\n\nTransforms the coordinates of the asymmetric unit cell to attend the translational symmetry of the cellulose phase (Iα, Iβ, II or III). This function returns the fractional coordinates of the unit cell.\n\n\n\n\n\n","category":"method"},{"location":"#CelluloseBuilder.z_expansion-Tuple{Vector{String}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, String}","page":"Home","title":"CelluloseBuilder.z_expansion","text":"z_expansion(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, zsize::Int64; phase=\"Iβ\")\n\nThis function is able to propagate the crystalline system across the z-axis.\n\n\n\n\n\n","category":"method"}]
}
