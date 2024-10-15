"""
    mergingPDBs(pdbfiles::Vector{String}; new_pdbfile=nothing)

    Aims to combine multipe PDBs on a new one, only changing the atom indexes.
    Tips: You can modify the segments, resids, chains, etc, then you can merge all the PDBs.
"""
function mergingPDBs(pdbfiles::Vector{String}; new_pdbfile=nothing)

    new_pdbfile = isnothing(new_pdbfile) ? tempname() * ".pdb" : new_pdbfile
    pdb = Base.open(new_pdbfile, "w")

    Base.write(pdb, "HEADER  1 CelluloseBuilder.jl - mergingPDBs function\n")

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

"""
    updatingPDB(psfname::String, pdbname::String, new_value::String; vmd_column="segid", vmd_selection="all", all=true, new_pdbname=nothing, vmd="vmd", vmd_debug=false)

    Aims to update a PSF and PDB file with a new value in a specific column. It is possible to adapt the selection such as on VMD using the flag `vmd_selection`
    and even export just the modification (`all = true`).

    # Arguments
    - psfname::String: the PSF file name.
    - pdbname::String: the PDB file name.
    - new_value::String: the new value to be updated.
    - vmd_column::String: the column to be updated (default: "segid").
    - vmd_selection::String: the selection to be updated (default: "all").
    - all::Bool: if you want to export all the PDB file or just the selection (default: true).
    - new_pdbname::String: the new PDB file name (default: nothing).
    - vmd::String: the VMD executable (default: "vmd").
    - vmd_debug::Bool: if you want to see the VMD command (default: false).

    # Examples
    ```julia
    updatingPDB("/tmp/jl_heF2BQrW5M.psf", "/tmp/jl_heF2BQrW5M.pdb", "BGC", vmd_column="resname", vmd_selection="resname BGLC")
    updatingPDB("/tmp/jl_heF2BQrW5M.psf", "/tmp/jl_heF2BQrW5M.pdb", "BGC", vmd_column="resname", vmd_selection="resname BGLC", all=false)
    ```
"""
function updatingPDB(psfname::String, pdbname::String, new_value::String; vmd_column="segid", vmd_selection="all", all=true, new_pdbname=nothing, vmd="vmd", vmd_debug=false)
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname
    new_psfname = replace(new_pdbname, ".pdb" => ".psf")

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"$vmd_selection\"]\n")
    Base.write(vmdinput, "\$sel set $vmd_column \"$new_value\"\n")
    Base.write(vmdinput, "\n")
    if vmd_selection != "all" && all
        Base.write(vmdinput, "[atomselect top \"all\"] writepsf $new_psfname\n")
        Base.write(vmdinput, "[atomselect top \"all\"] writepdb $new_pdbname\n")
    else
        Base.write(vmdinput, "\$sel writepsf $new_psfname\n")
        Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    end

    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    if vmd_debug
        return vmdoutput
    else
        return new_psfname, new_pdbname
    end

end

"""
    updatingPDB(pdbname::String, new_value::String; vmd_column="segid", vmd_selection="all", all=true, new_pdbname=nothing, vmd="vmd", vmd_debug=false)

    Aims to update a PSF and PDB file with a new value in a specific column. It is possible to adapt the selection such as on VMD using the flag `vmd_selection`
    and even export just the modification (`all = true`).

    # Arguments
    - pdbname::String: the PDB file name.
    - new_value::String: the new value to be updated.
    - vmd_column::String: the column to be updated (default: "segid").
    - vmd_selection::String: the selection to be updated (default: "all").
    - all::Bool: if you want to export all the PDB file or just the selection (default: true).
    - new_pdbname::String: the new PDB file name (default: nothing).
    - vmd::String: the VMD executable (default: "vmd").
    - vmd_debug::Bool: if you want to see the VMD command (default: false).
"""
function updatingPDB(pdbname::String, new_value::String; vmd_column="segid", vmd_selection="all", all=true, new_pdbname=nothing, vmd="vmd", vmd_debug=false)
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $pdbname\n")
    Base.write(vmdinput, "set sel [atomselect top \"$vmd_selection\"]\n")
    Base.write(vmdinput, "\$sel set $vmd_column \"$new_value\"\n")
    Base.write(vmdinput, "\n")
    if vmd_selection != "all" && all
        Base.write(vmdinput, "[atomselect top \"all\"] writepdb $new_pdbname\n")
    else
        Base.write(vmdinput, "\$sel writepdb $new_pdbname\n")
    end
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    if vmd_debug
        return vmdoutput
    else
        return new_pdbname
    end

end

"""
    fibril_segid(psfname::String, pdbname::String, old_starting::String, new_starting::String; new_pdbname=nothing, vmd="vmd", vmd_debug=false)

    Updates the segid of a fibril on a PSF and PDB file. This function is useful to build the cellulose macrofibril, avoinding replicas of segments
    on the same PSF. For example, with to microfibrils with the chain called M1 on segment, it is possible to change to A1 ... B1, and so on.
"""
function fibril_segid(psfname::String, pdbname::String, old_starting::String, new_starting::String; new_pdbname=nothing, vmd="vmd", vmd_debug=false)
  
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

    if vmd_debug
        return vmdoutput
    else
        return new_psfname, new_pdbname
    end

end


"""
    leftoverPDBs(pdbname::String; resname="TIP3", segname="WAT")

    Separate atoms from a PDB file based on their `resname`s and create new PDBs based on a raw `segname` + chain value.
    The main application is to correct the water PSF/PDB due to the lack of resnum (max 9999).
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

function pdb_replacement(pdbname::String; new_pdbname=nothing, pattern="ATOM  ", replacement="HETATM")
    
    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname

    Base.open(pdbname, "r") do infile
    
        Base.open(new_pdbname, "w") do outfile
            for line in eachline(infile)
                modified_line = replace(line, pattern => replacement)
                Base.write(outfile, "$modified_line\n")
            end
        end

    end

    return new_pdbname

end


