function cellulose(
                a::Int64, b::Int64, c::Int64;
                phase="Iβ", pbc=nothing, covalent=true, topology_file=generate_cellulose_topology(),
                vmd="vmd", vmdDebug=false
            )

    if a <= 1 || b <= 1 || c < 1
        error("The `a` ($a) and `b` ($b) dimensions must be ≥ 1, but `c` must be > 1.")
    end

    valid_pbc = isnothing(pbc) || in(lowercase(pbc), Set(["a", "b", "all"]))
    if !valid_pbc
        error("The periodic boundary conditions must be setted as `a`, `b`, or `all` (both directions). If you don't want to apply, just let `pbc=nothing`.")
    end
    
    println("""

    Building a cellulose crystal with (a, b, c) unit cell pattern.
    This crystal is built with the $phase phase and the periodic covalent bonding is setted as $covalent.
    """)
    xyzsizes, lattice, basisvectors = unitcell(a, b, c, phase)
    
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
    iii. screening the chains based on $(nfragments) fragments of the raw crystal, respecting the fibril xy-plane restricion for $phase phase.
    """)
    new_xyzname, chains = pbcXYZ(xyzname, nfragments, xyzsizes, phase=phase, vmd=vmd)

    println("""
    iv.  generating the PSF/PDB files:
         - writing the PDBs for each of those $(length(chains)) chains.")
         - cleaning each fragment PDB.")
         - using the CHARMM topology file to build the final PDB/PSF with the fragments
    """)
    pdbnames = fragPDBs(new_xyzname, fragments=chains, vmd=vmd)
    for pdb in pdbnames
        cleanPDB!(pdb, charmm_atomstyle)
    end

    monomers = if in(lowercase(phase), Set(["ib", "iβ", "ii", "ia", "iα"]))
        2 * xyzsizes[3]
    else
        xyzsizes[3]
    end

    vmdoutput = if lowercase(phase) == "ii"
        writePDB(
            monomers, pdbnames, pdbname=replace(new_xyzname, ".xyz" => ".pdb"),
            phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file,
            hasInversion=true
        )
    else
        writePDB(
            monomers, pdbnames, pdbname=replace(new_xyzname, ".xyz" => ".pdb"),
            phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file
        )
    end
    
    cleaning_tmpfiles(replace(new_xyzname, ".xyz" => ""))
    
    println("""

    It's done!
    """)

    if vmdDebug
        return vmdoutput
    end
end

function cellulose(
                chains::Int64, monomers::Int64;
                phase="Iα", layer="m", covalent=true, topology_file=generate_cellulose_topology(),
                vmd="vmd", vmdDebug=false
            )

    valid_fibril = (monomers >= 2) || (monomers%2 == 0)
    if !valid_fibril
        error("The number of cellobiose units must be an integer value equal or greater than 1. The actual value is $(Int64(monomers/2)).")
    end

    types = Set(["monolayer", "m", "center", "c", "origin", "o"])
    if !in(layer, types)
        error("The monolayer type must be `monolayer`, `center`, or `origin`.")
    end

    println("""

    Building a cellulose monolayer with $monomers cellobiose units and $chains polymeric chains.
    This layer is built with the $phase phase and the periodic covalent bonding is setted as $covalent.
    """)
    xyzsizes, lattice, basisvectors = unitcell(layer, chains, monomers, phase)
    
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
    iii. screening the chains based on $(nfragments) fragments of the raw crystal, respecting the fibril xy-plane restricion for $phase phase.
    """)
    new_xyzname, chains = pbcXYZ(xyzname, xyzsizes, phase=phase, structure=layer, vmd=vmd)

    println("""
    iv.  generating the PSF/PDB files:
         - writing the PDBs for each of those $(length(chains)) chains.")
         - cleaning each fragment PDB.")
         - using the CHARMM topology file to build the final PDB/PSF with the fragments
    """)
    pdbnames = fragPDBs(new_xyzname, fragments=chains, vmd=vmd)
    for pdb in pdbnames
        cleanPDB!(pdb, charmm_atomstyle)
    end

    vmdoutput = writePDB(
                    monomers, pdbnames, pdbname=replace(new_xyzname, ".xyz" => ".pdb"),
                    phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file
                )
    
    cleaning_tmpfiles(replace(new_xyzname, ".xyz" => ""))
    
    println("""

    It's done!
    """)

    if vmdDebug
        return vmdoutput
    end
end

"""
    cellulose(monomers::Int64; phase="Iβ", fibril=nothing, ...)

Builds a cellulose fibril with a given number of cellobiose units based on the number of `monomers`.
"""
function cellulose(
                monomers::Int64;
                phase="Iβ", fibril=nothing, covalent=true, topology_file=generate_cellulose_topology(),
                vmd="vmd", vmdDebug=false
            )

    valid_fibril = (monomers >= 2) || (monomers%2 == 0)
    if !valid_fibril
        error("The number of cellobiose units must be an integer value equal or greater than 1. The actual value is $(Int64(monomers/2)).")
    end

    println("""

    Building a cellulose fibril with $monomers cellobiose units.
    This fibril is built with the $phase phase and the periodic covalent bonding is setted as $covalent.
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
    iii. screening the chains based on $(nfragments) fragments of the raw crystal, respecting the fibril xy-plane restricion for $phase phase.
    """)
    new_xyzname, chains = pbcXYZ(xyzname, xyzsizes, phase=phase, structure="f", fibril=fibril, vmd=vmd)

    println("""
    iv.  generating the PSF/PDB files:
         - writing the PDBs for each of those $(length(chains)) chains.")
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

    vmdoutput = writePDB(
                    monomers, pdbnames, pdbname=replace(new_xyzname, ".xyz" => ".pdb"),
                    phase=phase, covalent=covalent, vmd=vmd, topology_file=topology_file
                )
    
    cleaning_tmpfiles(replace(new_xyzname, ".xyz" => ""))
    
    println("""

    It's done!
    """)

    if vmdDebug
        return vmdoutput
    end
end






