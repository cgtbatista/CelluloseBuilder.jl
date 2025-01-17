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
##    basisvectors = lattice2basis(lattice, phase)
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
##    xyzsize = PBC([xsize, ysize, zsize], phase, pbc=pbc)
##    println("       + appling transformations on fractional coordinates needed for the phase $phase.")
##    println("       + transforming the asymmetric unit to the cartesian coordinates for every [a,b,c] = [$xsize,$ysize,$zsize] Å.")
##    xinit, yinit, zinit = fractional2cartesian(xyzsize, phase)
##    println("       + atomic labels for $phase.")
##    atomsinit, atomstype = atomnames(phase, length(xinit))
##    println("")
##
##    println("   2 - Extending the cellulose modifications of the atoms:")
##    println("       + cleaning the coordinates and atomic labels for $phase.")
##    atomsclean, xclean, yclean, zclean = xy_pruning(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
##    println("       + expanding the z coordinates for $phase.")
##    atomsexpnd, xexpnd, yexpnd, zexpnd = z_expansion(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
##    println("       + picking the number of fragments of the basic structure.")
##    xyzfile, vmdoutput = writeXYZ(atomsexpnd, xexpnd, yexpnd, zexpnd)
##    n_fragments = picking_fragments(vmdoutput)
##    println("")
##
##    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
##    vmdxyz, frag_sel, frag_units, vmdoutput2 = pbcXYZ(n_fragments, xyzsize[1], xyzsize[2], phase=phase, pbc=pbc, xyzfile=xyzfile, vmd=vmd)
##    println("")
##
##    println("   4 - Generating the PSF/PDB files:")    
##    println("       + writing the PDBs for each of those $frag_units fragment units.")
##    pdb_basename = fragPDBs(vmdxyz, frag_sel, frag_units, vmd=vmd)[2]
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
##        vmdoutput3 = writePDB(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file, check_inversion=true)
##    else
##        vmdoutput3 = writePDB(monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
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
##    xyzsize, lattice = PBC(monolayer, nchains, ncellobiose, phase)
##    xsize, ysize, zsize = xyzsize[1], xyzsize[2], xyzsize[3]
##
##    ## DEALING W/ UNIT CELLS -----------------------------------------------------------------
##    println("   ... BASIS VECTORS")
##    basisvectors = lattice2basis(lattice, phase)
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
##    xinit, yinit, zinit = fractional2cartesian(xyzsize, phase)
##    println("       + atomic labels for $phase.")
##    atomsinit, atomstype = atomnames(phase, length(xinit))
##    println("")
##
##    println("   2 - Extending the cellulose modifications of the atoms:")
##    println("       + cleaning the coordinates and atomic labels for $phase.")
##    atomsclean, xclean, yclean, zclean = xy_pruning(atomsinit, xinit, yinit, zinit, xyzsize, phase=phase)
##    println("       + expanding the z coordinates for $phase.")
##    atomsexpnd, xexpnd, yexpnd, zexpnd = z_expansion(atomsclean, xclean, yclean, zclean, xyzsize[3], phase=phase)
##    println("       + picking the number of fragments of the basic structure.")
##    xyzfile, vmdoutput = writeXYZ(atomsexpnd, xexpnd, yexpnd, zexpnd)
##    n_fragments = picking_fragments(vmdoutput)
##    println("")
##
##    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
##    vmdxyz, frag_sel, frag_units, vmdoutput2 = pbcXYZ(monolayer, xyzsize, phase=phase, xyzfile=xyzfile, vmd=vmd)
##    println("")
##
##    println("   4 - Generating the PSF/PDB files:")    
##    println("       + writing the PDBs for each of those $frag_units fragment units.")
##    pdb_basename = fragPDBs(vmdxyz, frag_sel, frag_units, vmd=vmd)[2]
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
##    vmdoutput3 = writePDB(n_monomers, tmpfragments, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
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
function cellulosebuilder(monomers::Int64; phase="Iβ", fibril=nothing, covalent=true, vmd="vmd", topology_file=generate_cellulose_topology(), vmdDebug=false)

    valid_fibril = (monomers >= 2) || (monomers%2 == 0)
    if !valid_fibril
        error("The number of cellobiose units must be an integer value equal or greater than 1. The actual value is $(Int64(monomers/2)).")
    end

    println("""

    Building a cellulose fibril with $monomers cellobiose units.
    This fibrils are built with the $phase phase and the periodic covalent bonding is setted as $covalent.
    """)
    xyzsizes, lattice, basisvectors = unitcell(monomers, phase)

    println("""
    The unit cell has the lattice of ($(lattice[1]), $(lattice[2]), $(lattice[3])) and its basis vectors are:

        a = $(round.(basisvectors[1], digits=1)); b = $(round.(basisvectors[2], digits=1)); c = $(round.(basisvectors[3], digits=1))


    Building the cellulose structure file
    -------------------------------------


    i.   getting the initial unit cell coordinates and atomic labels:
         - imposing the default translational symmetry for pbc.
         - transforming the ASU to cartesian coordinates for every [xsize, ysize, zsize] = [$(xyzsizes[1]),$(xyzsizes[2]),$(xyzsizes[3])] Å.
         - atomic labels for $phase.
    """)
    xyz, charmm_atomstyle = crystal(xyzsizes, phase)

    println("""
    ii.  extending the cellulose modifications of the atoms:
         - cleaning the coordinates and atomic labels for $phase.
         - expanding the z coordinates for $phase.
         - picking the number of fragments of the basic structure.
    """)
    xyzname, nfragments = rawXYZ(xyz, xyzsizes, phase=phase)
    
    println("""
    iii. screening the $(nfragments) fragments of the $xyzname crystal that respect the fibril xy-plane restricion for $phase phase.
    """)
    new_xyzname, chains = pbcXYZ(xyzname, xyzsizes, phase=phase, fibril=fibril, vmd=vmd)

    println("""
    iv.  generating the PSF/PDB files:
         - writing the PDBs for each of those $(length(chains)) fragment units.")
         - cleaning each fragment PDB.")
         - using the CHARMM topology file to build the final PDB/PSF with the fragments
    """)
    pdbnames = fragPDBs(new_xyzname, fragments=chains, vmd=vmd)
    for pdb in pdbnames
        cleanPDB!(pdb, charmm_atomstyle)
    end

    #n = if in(lowercase(phase), Set(["ib", "iβ", "ii"]))
    #    2 * xyzsizes[3]
    #else
    #    xyzsizes[3]
    #end

    vmdoutput = writePDB(monomers, pdbnames, phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file)
    
    cleaning_tmpfiles("cellulose")
    
    println("""

    That's all! :D
    """)

    if vmdDebug
        return vmdoutput
    end
end






