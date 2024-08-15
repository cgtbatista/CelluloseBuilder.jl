
"""

    cellulosebuilder(a::Int64, b::Int64, c::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=DEFAULT_CARB_TOPOLOGY_FILE)
    cellulosebuilder(monolayer::String, units::Int64, ncellobiose::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=DEFAULT_CARB_TOPOLOGY_FILE)
    cellulosebuilder(type::String, ncellobiose::Int64; phase="I-BETA", pbc=:nothing, pcb_c=:true)

Cellulose-builder builds Cartesian coordinates files for cellulose crystalline domains and plant cell wall cellulose elementary fibrils in PDB format.

## Arguments

# Method I
- `a::Int64`, `b::Int64`, `c::Int64`: The dimensions of the unit cell in Angstroms. The `a` is the number of glucose units along x axis,
  `b` is the number of glucose units along y axis, and `c` is the number of glucose units along z axis.

# Method II
- `monolayer::String`: The type of. It could be `center`, `origin`, or `monolayer`.
- `units::Int64`: The number of cellulose units in the monolayer.
- `ncellobiose::Int64`: The number of cellobiose units along the cellulose unit.

# Method III
- `type::String`: The type of cellulose to be built. It could be `crystalline` or `fibril`.
- `ncellobiose::Int64`: The number of cellobiose units along the cellulose unit.

# Default arguments
- `phase`: The cellulose phase to be built. The default is `I-BETA`, but it could be.
- `pbc`: The periodic boundary conditions to be applied. It's could be around `A`, `B`, or both directions `ALL`.
- `pcb_c`: The covalent bond throught the periodic boundary conditions to be applied in the c axis. It will linking the first residue to the last one.

### Examples

```jldoctest

julia > 

```

"""

function cellulosebuilder(a::Int64, b::Int64, c::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=DEFAULT_CARB_TOPOLOGY_FILE)

    if a <= 1 || b <= 1 || c < 1
        error("The dimensions of the unit cell must be greater than 1. Only c can be 1.")
    end
    if !isnothing(pbc) && pbc != :A && pbc != :B && pbc != :ALL && pbc != :a && pbc != :b && pbc != :all && pbc != :All
        error("The periodic boundary conditions must be setted as A, B, or ALL (both directions). If you don't want to apply, just let it as nothing.")
    end
    
    ## HEADER   ------------------------------------------------------------------------------
    println("")
    println("GENERETING $phase CELLULOSE WITH ($a, $b, $c) LATTICE DIMENSIONS!")
    if covalent; println("COVALENT TURNNED ON -- CONSIDERING THE PERIODIC COVALENT BONDING ACROSS THE BOX BORDERS..."); end
    println("")
    println("")
    xsize, ysize, zsize = a, b, c
    if phase == "Iβ" || phase == "Ib" || phase == "II" || phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        lattice = [a, b, c]
    end
    if phase == "Iα" || phase == "Ia"
        lattice = [c, b, a]
    end
    

    ## DEALING W/ UNIT CELLS -----------------------------------------------------------------
    println("   ... BASIS VECTORS")
    basisvectors = gettingBasisVectors(lattice, phase)
    println("")
    println("           X -- $(round.(basisvectors[1], digits=1))")
    println("           Y -- $(round.(basisvectors[2], digits=1))")
    println("           Z -- $(round.(basisvectors[3], digits=1))")
    println("")

    ## EXPORTING THE CELLULOSE UNIT ---------------------------------------------------------
    println("   ... PREPARING THE PDB AND PSF FILES")
    println("")

    println("   1 - Getting the initial unit cell coordinates and atomic labels:")
    println("       + imposing translational symmetry for $pbc.")
    xyzsize = gettingPBC([xsize, ysize, zsize], phase, pbc=pbc)
    println("       + appling transformations on fractional coordinates needed for the phase $phase.")
    println("       + transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.")
    xinit, yinit, zinit = unitcell2cartesian(xyzsize, phase)
    println("       + atomic labels for $phase.")
    atomsinit, atomstype = atomsvecString(phase, length(xinit))
    println("")

    println("   2 - Extending the cellulose modifications of the atoms:")
    println("       + cleaning the coordinates and atomic labels for $phase.")
    atomsclean, xclean, yclean, zclean = _XY_trimming_coords(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
    println("       + expanding the z coordinates for $phase.")
    atomsexpnd, xexpnd, yexpnd, zexpnd = _Z_propagation_coords(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
    println("       + picking the number of fragments of the basic structure.")
    xyzfile, vmdoutput = _exporting_XYZfile(atomsexpnd, xexpnd, yexpnd, zexpnd)
    n_fragments = picking_fragments(vmdoutput)
    println("")

    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
    vmdxyz, frag_sel, frag_units, vmdoutput2 = transformingPBC(n_fragments, xyzsize[1], xyzsize[2], phase=phase, pbc=pbc, xyzfile=xyzfile, vmd=vmd)
    println("")

    println("   4 - Generating the PSF/PDB files:")    
    println("       + writing the PDBs for each of those $frag_units fragment units.")
    pdb_basename = _XYZfragments_2_PDB(vmdxyz, frag_sel, vmd=vmd)[2]
    println("       + cleaning each fragment PDB.")
    units = Base.split(frag_sel, " ");
    tmpfragments = String[]; tmpfile = tempname();
    for u in units
        pdbname = pdb_basename * "_" * u * ".pdb"
        new_pdbname = tmpfile * "_" * u * ".pdb"
        _cleaning_PDBfragment(atomstype, pdbname, new_pdbname)
        push!(tmpfragments, new_pdbname)
    end
    println("       + using the CHARMM topology file to build the final PDB/PSF with the fragments")
    if phase == "Iβ" || phase == "Ib" || phase == "II" || phase == "Iα" || phase == "Ia"
        monomers = 2*xyzsize[3]
    else monomers = xyzsize[3] end ### O III_I tá certo???? Contar o número de monômeros em cada cadeia depois...

    if phase == "II"
        vmdoutput3 = _exporting_PDBfile(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file, check_inversion=true)
    else
        vmdoutput3 = _exporting_PDBfile(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
    end
    
    cleaning_tmpfiles("cellulose")
    println("")
    println("   ... it is done!")

    return vmdoutput3

end


function cellulosebuilder(monolayer::String, units::Int64, ncellobiose::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=DEFAULT_CARB_TOPOLOGY_FILE)

    if monolayer != "monolayer" && monolayer != "center" && monolayer != "origin"
        error("The monolayer type must be `monolayer`, `center`, or `origin`.")
    end
    if units < 1 || ncellobiose < 1
        error("The number of units and cellobiose must be greater than 1.")
    end
    if !isnothing(pbc)
        pbc=nothing
    end
    
    ## HEADER   ------------------------------------------------------------------------------
    println("")
    println("GENERETING $phase CELLULOSE MONOLAYER ($monolayer)!")
    if covalent; println("COVALENT TURNNED ON -- CONSIDERING THE PERIODIC COVALENT BONDING ACROSS THE BOX BORDERS..."); end
    println("")
    println("")
    xyzsize, lattice = gettingPBC(monolayer, units, ncellobiose, phase)
    xsize, ysize, zsize = xyzsize[1], xyzsize[2], xyzsize[3]

    ## DEALING W/ UNIT CELLS -----------------------------------------------------------------
    println("   ... BASIS VECTORS")
    basisvectors = gettingBasisVectors(lattice, phase)
    println("")
    println("           X -- $(round.(basisvectors[1], digits=1))")
    println("           Y -- $(round.(basisvectors[2], digits=1))")
    println("           Z -- $(round.(basisvectors[3], digits=1))")
    println("")

    ## EXPORTING THE CELLULOSE UNIT ---------------------------------------------------------
    println("   ... PREPARING THE PDB AND PSF FILES")
    println("")

    println("   1 - Getting the initial unit cell coordinates and atomic labels:")
    println("       + imposing translational symmetry for $pbc.")
    println("       + transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.")
    xinit, yinit, zinit = unitcell2cartesian(xyzsize, phase)
    println("       + atomic labels for $phase.")
    atomsinit, atomstype = atomsvecString(phase, length(xinit))
    println("")

    println("   2 - Extending the cellulose modifications of the atoms:")
    println("       + cleaning the coordinates and atomic labels for $phase.")
    atomsclean, xclean, yclean, zclean = _XY_trimming_coords(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
    println("       + expanding the z coordinates for $phase.")
    atomsexpnd, xexpnd, yexpnd, zexpnd = _Z_propagation_coords(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
    println("       + picking the number of fragments of the basic structure.")
    xyzfile, vmdoutput = _exporting_XYZfile(atomsexpnd, xexpnd, yexpnd, zexpnd)
    n_fragments = picking_fragments(vmdoutput)
    println("")

    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
    vmdxyz, frag_sel, frag_units, vmdoutput2 = transformingPBC(monolayer, xyzsize, phase=phase, xyzfile=xyzfile, vmd=vmd)
    println("")

    println("   4 - Generating the PSF/PDB files:")    
    println("       + writing the PDBs for each of those $frag_units fragment units.")
    pdb_basename = _XYZfragments_2_PDB(vmdxyz, frag_sel, vmd=vmd)[2]
    println("       + cleaning each fragment PDB.")
    units = Base.split(frag_sel, " ");
    tmpfragments = String[]; tmpfile = tempname();
    for u in units
        pdbname = pdb_basename * "_" * u * ".pdb"
        new_pdbname = tmpfile * "_" * u * ".pdb"
        _cleaning_PDBfragment(atomstype, pdbname, new_pdbname)
        push!(tmpfragments, new_pdbname)
    end
    println("       + using the CHARMM topology file to build the final PDB/PSF with the fragments")
    if phase == "Iβ" || phase == "Ib" || phase == "II"
        monomers = 2*xyzsize[3]
    else monomers = xyzsize[3] end
    vmdoutput3 = _exporting_PDBfile(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
    
    cleaning_tmpfiles("cellulose")
    println("")
    println("   ... it is done!")

    return vmdoutput3

end









