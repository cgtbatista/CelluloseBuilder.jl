##function cellulosebuilder(a::Int64, b::Int64, c::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=generate_cellulose_topology())
##
##    if a <= 1 || b <= 1 || c < 1
##        error("The dimensions of the unit cell must be greater than 1. Only c can be 1.")
##    end
##    if !isnothing(pbc) && pbc != :A && pbc != :B && pbc != :ALL && pbc != :a && pbc != :b && pbc != :all && pbc != :All
##        error("The periodic boundary conditions must be setted as A, B, or ALL (both directions). If you don't want to apply, just let it as nothing.")
##    end
##    
##    ## HEADER   ------------------------------------------------------------------------------
##    println("")
##    println("GENERETING $phase CELLULOSE WITH ($a, $b, $c) LATTICE DIMENSIONS!")
##    if covalent; println("COVALENT TURNNED ON -- CONSIDERING THE PERIODIC COVALENT BONDING ACROSS THE BOX BORDERS..."); end
##    println("")
##    println("")
##    xsize, ysize, zsize = a, b, c
##    if phase == "Iβ" || phase == "Ib" || phase == "II" || phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
##        lattice = [a, b, c]
##    end
##    if phase == "Iα" || phase == "Ia"
##        lattice = [c, b, a]
##    end
##    
##
##    ## DEALING W/ UNIT CELLS -----------------------------------------------------------------
##    println("   ... BASIS VECTORS")
##    basisvectors = gettingBasisVectors(lattice, phase)
##    println("")
##    println("           X -- $(round.(basisvectors[1], digits=1))")
##    println("           Y -- $(round.(basisvectors[2], digits=1))")
##    println("           Z -- $(round.(basisvectors[3], digits=1))")
##    println("")
##
##    ## EXPORTING THE CELLULOSE UNIT ---------------------------------------------------------
##    println("   ... PREPARING THE PDB AND PSF FILES")
##    println("")
##
##    println("   1 - Getting the initial unit cell coordinates and atomic labels:")
##    println("       + imposing translational symmetry for $pbc.")
##    xyzsize = gettingPBC([xsize, ysize, zsize], phase, pbc=pbc)
##    println("       + appling transformations on fractional coordinates needed for the phase $phase.")
##    println("       + transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.")
##    xinit, yinit, zinit = unitcell2cartesian(xyzsize, phase)
##    println("       + atomic labels for $phase.")
##    atomsinit, atomstype = atomnames(phase, length(xinit))
##    println("")
##
##    println("   2 - Extending the cellulose modifications of the atoms:")
##    println("       + cleaning the coordinates and atomic labels for $phase.")
##    atomsclean, xclean, yclean, zclean = _trimming_xy(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
##    println("       + expanding the z coordinates for $phase.")
##    atomsexpnd, xexpnd, yexpnd, zexpnd = _expanding_z(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
##    println("       + picking the number of fragments of the basic structure.")
##    xyzfile, vmdoutput = writeXYZ(atomsexpnd, xexpnd, yexpnd, zexpnd)
##    n_fragments = picking_fragments(vmdoutput)
##    println("")
##
##    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
##    vmdxyz, frag_sel, frag_units, vmdoutput2 = transformingPBC(n_fragments, xyzsize[1], xyzsize[2], phase=phase, pbc=pbc, xyzfile=xyzfile, vmd=vmd)
##    println("")
##
##    println("   4 - Generating the PSF/PDB files:")    
##    println("       + writing the PDBs for each of those $frag_units fragment units.")
##    pdb_basename = xyz2pdb(vmdxyz, frag_sel, frag_units, vmd=vmd)[2]
##    println("       + cleaning each fragment PDB.")
##    units = Base.split(frag_sel, " ");
##    tmpfragments = String[]; tmpfile = tempname();
##    for u in units
##        pdbname = pdb_basename * "_" * u * ".pdb"
##        new_pdbname = tmpfile * "_" * u * ".pdb"
##        cleanPDB(atomstype, pdbname, new_pdbname)
##        push!(tmpfragments, new_pdbname)
##    end
##    println("       + using the CHARMM topology file to build the final PDB/PSF with the fragments")
##    if phase == "Iβ" || phase == "Ib" || phase == "II" || phase == "Iα" || phase == "Ia"
##        monomers = 2*xyzsize[3]
##    else monomers = xyzsize[3] end ### O III_I tá certo???? Contar o número de monômeros em cada cadeia depois...
##
##    if phase == "II"
##        vmdoutput3 = _exporting_PDBfile(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file, check_inversion=true)
##    else
##        vmdoutput3 = _exporting_PDBfile(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
##    end
##    
##    cleaning_tmpfiles("cellulose")
##    println("")
##    println("   ... it is done!")
##
##    return vmdoutput3
##
##end
##
##function cellulosebuilder(monolayer::String, nchains::Int64, monomers::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd", topology_file=generate_cellulose_topology())
##
##    ncellobiose = Int64(monomers/2)
##
##    if monomers < 2
##        error("The number of cellobiose units must be equal or greater than 1. The actual value is $(ncellobiose).")
##    end ## checking the number of cellobiose units
##
##    if monomers%2 != 0
##        error("The number of monomer units must be an even number. The actual value is $(monomers%2).")
##    end ## checking for a even number of cellobiose units
##
##    types = Set(["monolayer", "center", "origin"])
##
##    if !in(monolayer, types)
##        error("The monolayer type must be `monolayer`, `center`, or `origin`.")
##    end
##
##    pbc = nothing
##    
##    ## HEADER   ------------------------------------------------------------------------------
##    println("")
##    println("GENERETING $phase CELLULOSE MONOLAYER ($monolayer)!")
##    if covalent; println("COVALENT TURNNED ON -- CONSIDERING THE PERIODIC COVALENT BONDING ACROSS THE BOX BORDERS..."); end
##    println("")
##    println("")
##    xyzsize, lattice = gettingPBC(monolayer, nchains, ncellobiose, phase)
##    xsize, ysize, zsize = xyzsize[1], xyzsize[2], xyzsize[3]
##
##    ## DEALING W/ UNIT CELLS -----------------------------------------------------------------
##    println("   ... BASIS VECTORS")
##    basisvectors = gettingBasisVectors(lattice, phase)
##    println("")
##    println("           X -- $(round.(basisvectors[1], digits=1))")
##    println("           Y -- $(round.(basisvectors[2], digits=1))")
##    println("           Z -- $(round.(basisvectors[3], digits=1))")
##    println("")
##
##    ## EXPORTING THE CELLULOSE UNIT ---------------------------------------------------------
##    println("   ... PREPARING THE PDB AND PSF FILES")
##    println("")
##
##    println("   1 - Getting the initial unit cell coordinates and atomic labels:")
##    println("       + imposing translational symmetry for $pbc.")
##    println("       + transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.")
##    xinit, yinit, zinit = unitcell2cartesian(xyzsize, phase)
##    println("       + atomic labels for $phase.")
##    atomsinit, atomstype = atomnames(phase, length(xinit))
##    println("")
##
##    println("   2 - Extending the cellulose modifications of the atoms:")
##    println("       + cleaning the coordinates and atomic labels for $phase.")
##    atomsclean, xclean, yclean, zclean = _trimming_xy(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
##    println("       + expanding the z coordinates for $phase.")
##    atomsexpnd, xexpnd, yexpnd, zexpnd = _expanding_z(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
##    println("       + picking the number of fragments of the basic structure.")
##    xyzfile, vmdoutput = writeXYZ(atomsexpnd, xexpnd, yexpnd, zexpnd)
##    n_fragments = picking_fragments(vmdoutput)
##    println("")
##
##    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
##    vmdxyz, frag_sel, frag_units, vmdoutput2 = transformingPBC(monolayer, xyzsize, phase=phase, xyzfile=xyzfile, vmd=vmd)
##    println("")
##
##    println("   4 - Generating the PSF/PDB files:")    
##    println("       + writing the PDBs for each of those $frag_units fragment units.")
##    pdb_basename = xyz2pdb(vmdxyz, frag_sel, frag_units, vmd=vmd)[2]
##    println("       + cleaning each fragment PDB.")
##    units = Base.split(frag_sel, " ");
##    tmpfragments = String[]; tmpfile = tempname();
##    for u in units
##        pdbname = pdb_basename * "_" * u * ".pdb"
##        new_pdbname = tmpfile * "_" * u * ".pdb"
##        cleanPDB(atomstype, pdbname, new_pdbname)
##        push!(tmpfragments, new_pdbname)
##    end
##    println("       + using the CHARMM topology file to build the final PDB/PSF with the fragments")
##    if phase == "Iβ" || phase == "Ib" || phase == "II"
##        n_monomers = 2*xyzsize[3]
##    else n_monomers = xyzsize[3] end
##    vmdoutput3 = _exporting_PDBfile(n_monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
##    
##    cleaning_tmpfiles("cellulose")
##    println("")
##    println("   ... it is done!")
##
##    return vmdoutput3
##
##end

"""
    cellulosebuilder(monomers::Int64; phase="Iβ", fibril=nothing, ...)

Builds a cellulose fibril with a given number of cellobiose units based on the number of `monomers`.
"""
function cellulosebuilder(monomers::Int64; phase="Iβ", fibril=nothing, covalent=true, vmd="vmd", topology_file=generate_cellulose_topology())

    ncellobiose = Int64(monomers/2)

    if monomers < 2
        error("The number of cellobiose units must be equal or greater than 1. The actual value is $(ncellobiose).")
    end ## checking the number of cellobiose units

    if monomers%2 != 0
        error("The number of monomer units must be an even number. The actual value is $(monomers%2).")
    end ## checking for a even number of cellobiose units

    println("""

    Building a cellulose fibril with $monomers cellobiose units.
    This fibrils are built with the $phase phase and the periodic covalent bonding is setted as $covalent.
    """)
    xyzsize, lattice = gettingPBC(ncellobiose, phase)
    basisvectors = gettingBasisVectors(lattice, phase)

    println("""
    The (a, b, c) unit cell dimensions are ($(xsize = xyzsize[1]), $(ysize = xyzsize[2]), $(zsize = xyzsize[3])) and the basis vectors are:

        x = $(round.(basisvectors[1], digits=1)); y = $(round.(basisvectors[2], digits=1)); z = $(round.(basisvectors[3], digits=1))


    Building the cellulose structure file
    -------------------------------------
    """)

    pbc = nothing

    println("""
    i.   getting the initial unit cell coordinates and atomic labels:
         - imposing translational symmetry for pbc = $pbc.
         - transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.
         - atomic labels for $phase.
    """)
    x, y, z = unitcell2cartesian(xyzsize, phase)
    atoms, atomstype = atomnames(phase, nblocks=length(x))
    crystal = XYZs(atoms, x, y, z)

    println("""
    ii.  extending the cellulose modifications of the atoms:
         - cleaning the coordinates and atomic labels for $phase.
         - expanding the z coordinates for $phase.
         - picking the number of fragments of the basic structure.
    """)
    xyzfile, nfragments = rawXYZ(crystal, xyzsize, phase=phase)
    
    println("""
    iii. applying the periodic boundary conditions (PBC) on the $(nfragments) fragments: $(pbc)...
    """)
    vmdxyz, frag_sel, frag_units = transformingPBC("fibril", xyzsize, phase=phase, fibril=fibril, xyzfile=xyzfile, vmd=vmd)

    println("""
    iv.  generating the PSF/PDB files:
         - writing the PDBs for each of those $frag_units fragment units.")
         - cleaning each fragment PDB.")
         - using the CHARMM topology file to build the final PDB/PSF with the fragments
    """)
    pdb_basename = xyz2pdb(vmdxyz, frag_sel, frag_units, vmd=vmd)
    units = Base.split(frag_sel, " ");
    tmpfragments = String[]; tmpfile = tempname();
    for u in units
        pdbname = pdb_basename * "_" * u * ".pdb"
        new_pdbname = tmpfile * "_" * u * ".pdb"
        cleanPDB(atomstype, pdbname, new_pdbfile=new_pdbname)
        push!(tmpfragments, new_pdbname)
    end
    if phase == "Iβ" || phase == "Ib" || phase == "II"
        n_monomers = 2*xyzsize[3]
    else n_monomers = xyzsize[3] end
    vmdoutput3 = _exporting_PDBfile(n_monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
    
    cleaning_tmpfiles("cellulose")
    
    println("""

    That's all! :D
    """)

    return vmdoutput3

end






