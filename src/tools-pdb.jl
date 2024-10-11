"""
    mergingPDBs(pdbfiles::Vector{String}; new_pdbfile=nothing)

    Aims to combine multipe PDBs on a new one, only changing the atom indexes.
"""
function mergingPDBs(pdbfiles::Vector{String}; new_pdbfile=nothing)

    new_pdbfile = isnothing(new_pdbfile) ? tempname() * ".pdb" : new_pdbfile
    pdb = Base.open(new_pdbfile, "w")

    Base.write(pdb, "REMARK  1 CelluloseBuilder.jl - mergingPDBs function\n")

    atom_index = 0

    for pdbfile in pdbfiles

        # editing a temporary PDB

        tmpPDB = PDBTools.readPDB(pdbfile)

        for at in eachindex(tmpPDB)
            atom_index += 1
            tmpPDB[at].index = atom_index
            tmpPDB[at].index_pdb = atom_index
        end

        PDBTools.writePDB(tmpPDB, pdbfile)

        # writing some lines of the new PDB

        lines = readlines(pdbfile)
        content = lines[2:end-1]
        for line in content
            Base.write(pdb, "$line\n")
        end

    end

    Base.write(pdb, "END\n")

    Base.close(pdb)

    return new_pdbfile
    
end


##function new_mergingPDBs(pdbfiles::Vector{String}; filename=nothing)
##    
##    filename = isnothing(filename) ? tempname() * ".pdb" : filename
##
##    pdb = open(expanduser(filename), "w")
##
##    println(pdb, "CelluloseBuilder.jl - mergingPDBs function")
##    atom_index = 0
##    # atoms::AbstractVector{Atom}
##
##    for pdbfile in pdbfiles
##
##        atoms = PDBTools.readPDB(pdbfile)
##
##        for a in eachindex(atoms)
##            atom_index += 1
##            atoms[a].index = atom_index
##            atoms[a].index_pdb = atom_index
##            println(pdb, PDBTools.write_atom(atoms[a]))
##        end
##
##    end
##
##    println(pdb, "END")
##
##    close(pdb)
##
##    return filename
##end


function updating_resid(pdbname::String, resid::Int64; new_pdbname=nothing)
  
    pdb = PDBTools.readPDB(pdbname)  

    for at in pdb
        at.resnum = resid
    end

    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname
    PDBTools.writePDB(pdb, new_pdbname)

    return new_pdbname

end

function updating_segid(pdbname::String, segid::String; new_pdbname=nothing, vmd="vmd")
  
    if isnothing(new_pdbname)
        new_pdbname = tempname() * ".pdb"
    end

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"all\"]\n")
    Base.write(vmdinput, "\$sel set segid \"$segid\"\n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return vmdoutput

end

function updating_segid(psfname::String, pdbname::String, segid::String; new_name=nothing, vmd="vmd")
  
    if isnothing(new_name)
        new_psfname = tempname() * ".psf"
        new_pdbname = tempname() * ".pdb"
    else
        new_psfname = new_name * ".psf"
        new_pdbname = new_name * ".pdb"
    end

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"all\"]\n")
    Base.write(vmdinput, "\$sel set segid \"$segid\"\n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "\$sel writepsf $new_psfname\n")
    Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return vmdoutput

end

function updatingPDB(psfname::String, pdbname::String, new_value::String; vmd_column="segid", new_pdbname=nothing, vmd="vmd")
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname
    new_psfname = replace(new_pdbname, ".pdb" => ".psf")

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"all\"]\n")
    Base.write(vmdinput, "\$sel set $vmd_column \"$new_value\"\n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "\$sel writepsf $new_psfname\n")
    Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return new_psfname, new_pdbname, vmdoutput

end

function updatingPDB(pdbname::String, new_value::String; vmd_column="segid", new_pdbname=nothing, vmd="vmd")
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"all\"]\n")
    Base.write(vmdinput, "\$sel set $vmd_column \"$new_value\"\n")
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return new_pdbname, vmdoutput

end

function fibril_segid(psfname::String, pdbname::String, old_starting::String, new_starting::String; new_pdbname=nothing, vmd="vmd")
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname
    new_psfname = replace(new_pdbname, ".pdb" => ".psf")

    pdb_segnames = PDBTools.segname.(PDBTools.readPDB(pdbname))
    new_segnames = replace.(pdb_segnames, old_starting => new_starting)

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n")
    Base.write(vmdinput, "\n")
    
    for ith_segids in eachindex(new_segnames)
        Base.write(vmdinput, "[atomselect top \"index $(ith_segids-1)\"] set segid \"$(new_segnames[ith_segids])\"\n")
    end
    
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepsf $new_psfname\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return new_psfname, new_pdbname, vmdoutput

end


"""
"""
function leftoverPDBs(pdbname::String; resname="TIP3", segname="WAT")
    
    pdbfiles = []
    
    if typeof(resname) != Vector{String}
        resname = [resname]
    end

    if !isnothing(segname) && (typeof(segname) != Vector{String})
        segname = [segname]
    end

    if !isnothing(segname) && (length(resname) != length(segname))
        throw(ArgumentError("The length of the segname must be the same as the resname..."))
    end

    for ith_resname in eachindex(resname)

        pdb = PDBTools.readPDB(pdbname,
                    only = (atom -> atom.resname == resname[ith_resname])
                )
        
        chains = unique(PDBTools.chain.(pdb))

        for chain in chains
            chain_pdb = PDBTools.select(pdb,
                    by = (atom -> atom.chain == chain)
                ); chain_pdbname = tempname() * ".pdb"

            PDBTools.writePDB(chain_pdb, chain_pdbname)

            if isnothing(segname)
                new_segment = String("$(resname[ith_resname][1:3])$chain")
            else
                new_segment = String("$(segname[ith_resname])$chain")
            end

            push!(pdbfiles, (new_segment, chain_pdbname))
        end

    end

    return pdbfiles

end

function dummy_chains()
    return "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"
end

