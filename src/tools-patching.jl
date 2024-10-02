function patching(
    chain_pdbname::String, resid::Vector{Int64}, chain_name::String, segid::String;
    decoration="PETN", new_segid=nothing, new_resid=nothing, new_chain_name=nothing,
    filename=nothing, patches_pdbname=nothing,
    vmd="vmd", topology=DEFAULT_CARB_TOPOLOGY_FILE)

    if isnothing(filename)
        filename = tempname()
        psfname = filename * ".psf"
        pdbname = filename * ".pdb"
    end

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "package require psfgen\n") 
    if typeof(topology) == Vector{String}
        for top in topology
            Base.write(vmdinput, "topology $top\n")
        end
    else
        Base.write(vmdinput, "topology $topology\n")
    end
    Base.write(vmdinput, "\n")

    Base.write(vmdinput, "readpsf $(replace(chain_pdbname, ".pdb" => ".psf"))\n")
    Base.write(vmdinput, "coordpdb $chain_pdbname\n\n")

    patchings = Vector{Int64}[]; pdb_decorations = String[]
    pdb = PDBTools.readPDB(chain_pdbname)
    ith_patch = maximum(PDBTools.resnum.(pdb))
    for i in resid
        ith_patch += ith_patch
        new_patch = matching_residue(chain_pdbname, i, chain_name, segid, decoration=decoration, new_resid=ith_patch, new_chain=new_chain_name, new_segid=new_segid)
        push!(patchings, [i, ith_patch])
        push!(pdb_decorations, new_patch)
    end

    patches_pdb = mergingPDBs(pdb_decorations, new_pdbfile=patches_pdbname)

    Base.write(vmdinput, "segment $new_segid { pdb $patches_pdb }\n")

    for patch in patchings
        Base.write(vmdinput, "patch PCEL $segid:$(patch[1]) $new_segid:$(patch[2])\n")
    end
    Base.write(vmdinput, "\n")

    Base.write(vmdinput, "regenerate angles dihedrals\n")  
    Base.write(vmdinput, "coordpdb $patches_pdb $new_segid\n")
    Base.write(vmdinput, "guesscoord\n\n")
    Base.write(vmdinput, "writepsf $psfname\n")
    Base.write(vmdinput, "writepdb $pdbname\n\n")

    Base.write(vmdinput, "mol new     $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n\n")
    Base.write(vmdinput, "set sel [atomselect top \"all\"]\n")
    Base.write(vmdinput, "\$sel set segid \"$segid\"\n\n")
    Base.write(vmdinput, "\$sel writepsf $psfname\n\n")
    Base.write(vmdinput, "\$sel writepdb $pdbname\n\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return vmdoutput

end


function matching_residue(
        main_pdbname::String, resid::Int64, chain::String, segid::String;
        decoration="PETN", new_resid=nothing, new_chain=nothing, new_segid=nothing, decoration_pdbname=nothing, new_pdbname=nothing
    )
  
    decoration_pdbname = isnothing(decoration_pdbname) ? getPDB(decoration) : decoration_pdbname
    tmpPDB = PDBTools.readPDB(decoration_pdbname)

    idxs, ref0_coord, ref1_coord, ini0_coord, ini1_coord = decoration_library(main_pdbname, decoration_pdbname, decoration, resid=resid, segid=segid)

    # Rotating the residue coordinates following the β-Glc residue orientation...
    v_residue = ini1_coord - ini0_coord
    v_target = ref0_coord - ref1_coord
    rotated_coords = rotate_residue(
            PDBTools.coor.(tmpPDB),
            ref1_coord,
            R_matrix(v_target, v_residue)
        )
    # Translating the residue coordinates to the specify position...
    v_trans = ref0_coord - rotated_coords[idxs[3]]
    translated_coords = translate_residue(
            rotated_coords, v_trans
        )

    for at in eachindex(translated_coords)
        tmpPDB[at].x = translated_coords[at][1]
        tmpPDB[at].y = translated_coords[at][2]
        tmpPDB[at].z = translated_coords[at][3]
        tmpPDB[at].chain = ifelse(isnothing(new_chain), chain, new_chain)
        tmpPDB[at].resnum = ifelse(isnothing(new_resid), resid, new_resid)
    end

    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname
    PDBTools.writePDB(tmpPDB, new_pdbname)

    vmdoutput = updating_segid(
            new_pdbname,
            ifelse(isnothing(new_segid), segid, new_segid),
            new_pdbname=new_pdbname
        ) ## não sei o porquê, mas não consigo salvar o nome do segmento no PDBTools: `res[at].segname = ifelse(isnothing(new_segid), segid, new_segid)`

    return new_pdbname

end






