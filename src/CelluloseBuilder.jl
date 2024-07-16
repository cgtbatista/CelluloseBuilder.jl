module CelluloseBuilder

# Write your package code here.

"""

    cellulosebuilder(a::Int64, b::Int64, c::Int64; phase="I-BETA", pbc=:nothing, pcb_c=:true)
    cellulosebuilder(monolayer::String, units::Int64, ncellobiose::Int64; phase="I-BETA", pbc=:nothing, pcb_c=:true)
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

```julia-repl

julia > 

```

"""

function cellulosebuilder(a::Int64, b::Int64, c::Int64; phase="Iβ", pbc=nothing, covalent=true, vmd="vmd")

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
    lattice = [a, b, c]

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
    ## PHASES CAN DIFFER...
    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
        xinit, yinit, zinit = unitcell2cartesian(xyzsize[1:2], phase)
    end
    println("       + atomic labels for $phase.")
    atomsinit, atomstype = atomsvecString(phase, xyzsize[1], xyzsize[2])
    println("")

    println("   2 - Extending the cellulose modifications of the atoms:")
    println("       + cleaning the coordinates and atomic labels for $phase.")
    atomsclean, xclean, yclean, zclean = XY_coord_trimming(atomsinit, xinit, yinit, zinit, xyzsize, phase)
    println("       + expanding the z coordinates for $phase.")
    atomsexpnd, xexpnd, yexpnd, zexpnd = Z_propagation(atomsclean, xclean, yclean, zclean, xyzsize[3], phase)
    println("       + picking the number of fragments of the basic structure.")
    xyzfile, vmdoutput = _exporting_XYZfile(atomsexpnd, xexpnd, yexpnd, zexpnd)
    n_fragments = picking_fragments(vmdoutput)
    println("")

    println("   3 - Periodic boundary conditions (PBC) on the $n_fragments fragments: $(pbc)...")
    vmdxyz, frag_sel, frag_units = PBCtools(phase, pbc, n_fragments, xyzsize[1], xyzsize[2], xyzfile, vmd)
    println("")

    println("   4 - Generating the PSF/PDB files:")    
    println("       + writing the PDBs for each of those $frag_units fragment units.")
    pdb_basename = _XYZfragments_2_PDB(vmdxyz, frag_sel, vmd=vmd)[1]
    println("       + cleaning each fragment PDB.")
    units = Base.split(frag_sel, " ");
    tmpfragments = String[];
    for u in units
        pdbname = pdb_basename * "_" * u * ".pdb"
        new_pdbname = "/tmp/tmp_" * u * ".pdb"
        _PDBfragment_cleaning(atomstype, pdbname, new_pdbname)
        push!(tmpfragments, new_pdbname)
    end
    println("       + using the CHARMM topology file to build the final PDB/PSF with the fragments")
    #topology_file = "/home/geckodo/Documents/dupree/fibril/cellulose-builder/top_all36_carb.rtf"
    topology_file = "/home/user/Documents/phd/dupree/fibril/cellulose-builder/top_all36_carb.rtf"
    vmdoutput2 = _exporting_PDBfile(phase, 2*xyzsize[3], tmpfragments, topology_file, covalent=covalent, vmd=vmd)
    
    cleaning_tmpfiles()
    println("")
    println("   ... it is done!")

    return vmdoutput2
end

function _exporting_PDBfile(phase, max_monomer, tmp, topology_file; covalent=true, vmd="vmd")

    fraglist = collect(1:1:length(tmp))
    monomer = 0; prev_monomer = 0;

    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
            
        vmdinput = Base.open("/tmp/patching.tcl", "w")
        
        Base.write(vmdinput, "package require psfgen \n")
        Base.write(vmdinput, "topology $(topology_file) \n")
        Base.write(vmdinput, " \n")
        for frag in fraglist
            segname = "M$frag"
            Base.write(vmdinput, "segment $segname { \n")
            Base.write(vmdinput, "    pdb $(tmp[frag]) \n")
            Base.write(vmdinput, "} \n")
            monomer = max_monomer
            while monomer > 1
                prev_monomer = monomer - 1
                Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$prev_monomer \n")
                monomer -= 1
                if monomer == 1 && covalent == true
                    Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$max_monomer \n")
                end
            end
        end
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "regenerate angles dihedrals \n")
        Base.write(vmdinput, " \n")
        for frag in fraglist
            segname = "M$frag"
            Base.write(vmdinput, "coordpdb $(tmp[frag]) $segname \n")
        end
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "guesscoord \n")
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "writepsf /tmp/cellulose.psf \n")
        Base.write(vmdinput, "writepdb /tmp/cellulose.pdb \n")
        Base.close(vmdinput)
        vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e /tmp/patching.tcl`, String), "\n")
    end

    return vmdoutput
end

function _PDBfragment_cleaning(atoms::Vector{String}, pdbfile::String, new_pdbfile::String)

    ith_atom = 1; resid = 1;

    pdbdata = Base.split(Base.read(pdbfile, String), "\n")
    pdb = open(new_pdbfile, "w")

    Base.write(pdb, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")

    for irow in pdbdata
        if !occursin("ATOM  ", irow); continue; end            
        
        new_irow = irow
        atomname = atoms[ith_atom]
        
        # setting hetatoms
        new_irow = replace(new_irow, "ATOM  " => "HETATM")
        # editing the atom names
        _atomname_length = length(atomname); _pdbcolumn_name = SubString(irow, 14 : 14+_atomname_length-1)
        new_irow = replace(new_irow, _pdbcolumn_name => atomname, count=1)
        
        # editing the residue names
        _pdbcolumn_resname = SubString(new_irow, findfirst("    X ", new_irow))
        new_irow = replace(new_irow, _pdbcolumn_resname => "BGLC  ")
        # editing residue sequence numbers
        _pdbcolumn_resnum = SubString(new_irow, 26-length("$resid")+1 : 26)
        new_irow = replace(new_irow, _pdbcolumn_resnum*"    " => "$resid"*"    ") ## This dummy space is escential to not replace all numbers (x coords begin on 31 string position)
        # editing the segname
        _pdbcolumn_segname = SubString(new_irow, findfirst("  0.00  0.00           ", irow))
        new_irow = replace(new_irow, _pdbcolumn_segname => "  1.00  0.00           ")
        
        # writting the new configuration
        Base.write(pdb, "$new_irow\n")
        
        ith_atom += 1
        if ith_atom > length(atoms)
            ith_atom = 1
            resid += 1
        end
    end
    Base.write(pdb, "END")
    
    Base.close(pdb)

    return nothing

end

end
