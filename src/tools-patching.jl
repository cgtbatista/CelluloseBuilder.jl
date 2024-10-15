"""
   patching(
        reference::String, resid::Vector{Int64}, segid::String;
        decoration_type="ENP", new_segid="TMP", new_pdbname=nothing, patches_pdbname=nothing, vmd="vmd", vmd_debug=false,
        topology=[ generate_cellulose_topology(), generate_petn_topology() ]
    )

    This function patches a reference PDB file with decoration residues. The decoration residues are rotated and translated to match the reference residues.
    
    # Examples
    ```julia
    patching("cellulose.pdb", [1, 2, 3], "M1", decoration_type="ENP", new_segid="TMP", new_pdbname="cellulose_enp.pdb")
    ```
"""
function patching(
        reference::String, resid::Vector{Int64}, segid::String;
        decoration_type="ENP", new_segid="TMP", new_pdbname=nothing, patches_pdbname=nothing, vmd="vmd", vmd_debug=false,
        topology=[ generate_cellulose_topology(), generate_petn_topology() ]
    )

    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname
    new_psfname = replace(new_pdbname, ".pdb" => ".psf")

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "package require psfgen\n")
    if typeof(topology) == Vector{String}
        for topfile in topology
            Base.write(vmdinput, "topology $topfile\n")
        end
    else
        Base.write(vmdinput, "topology $topology\n")
    end
    Base.write(vmdinput, "\n")

    Base.write(vmdinput, "readpsf $(replace(reference, ".pdb" => ".psf"))\n")
    Base.write(vmdinput, "coordpdb $reference\n\n")

    patch_pairlist = Vector{Int64}[]; decorationPDBs = String[]

    last_reference_resid = maximum(
            PDBTools.resnum.(PDBTools.readPDB(reference))
        )

    for i in eachindex(resid)
        ith_chain_resid = resid[i]
        ith_decoration_resid = last_reference_resid + i
        new_patch = matching_residue(reference, ith_chain_resid, segid, decoration_type=decoration_type, new_resid=ith_decoration_resid)
        push!(patch_pairlist, [ith_chain_resid, ith_decoration_resid])
        push!(decorationPDBs, new_patch)
    end

    decorations_pdb = mergingPDBs(decorationPDBs, new_pdbfile=patches_pdbname)

    Base.write(vmdinput, "segment $new_segid { pdb $decorations_pdb }\n")

    for pair in patch_pairlist
        Base.write(vmdinput, "patch PCEL $segid:$(pair[1]) $new_segid:$(pair[2])\n")
    end
    Base.write(vmdinput, "\n")

    Base.write(vmdinput, "regenerate angles dihedrals\n")  
    Base.write(vmdinput, "coordpdb $decorations_pdb $new_segid\n")
    Base.write(vmdinput, "guesscoord\n\n")
    Base.write(vmdinput, "writepsf $new_psfname\n")
    Base.write(vmdinput, "writepdb $new_pdbname\n\n")

    ## Updating the decoration segname (is good that it have the same segname)
    Base.write(vmdinput, "mol new     $new_psfname\n")
    Base.write(vmdinput, "mol addfile $new_pdbname\n\n")
    Base.write(vmdinput, "set sel [atomselect top \"resname $decoration_type\"]\n")
    Base.write(vmdinput, "\$sel set segid \"$segid\"\n\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepsf $new_psfname\n\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepdb $new_pdbname\n\n")
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
   matching_residue(reference::String, resid::Int64, segid::String; decoration=nothing, decoration_type="ENP", new_resid=nothing, new_segid=nothing, new_pdbname=nothing)

    This function matches a residue from a reference PDB file to a decoration residue. The decoration residue is rotated and translated to match the reference residue.

    # Arguments
    - `reference::String`: The name of the reference PDB file.
    - `resid::Int64`: The residue number of the reference PDB file.
    - `segid::String`: The segment name of the reference PDB file.
    - `decoration::String`: The name of the decoration residue. Currently, only `ENP`.
    - `decoration_type::String`: The type of decoration residue. The default is `ENP`.
    - `new_resid::Int64`: The new residue number of the decoration residue. The default is `nothing`.
    - `new_segid::String`: The new segment name of the decoration residue. The default is `nothing`.
    - `new_pdbname::String`: The name of the new PDB file. The default is a temporary file.

    # Examples
    ```julia
    matching_residue("cellulose.pdb", 1, "M1", decoration="ENP", new_resid=2, new_segid="M2", new_pdbname="cellulose_enp.pdb")
    ```
"""
function matching_residue(
        reference::String, resid::Int64, segid::String;
        decoration=nothing, decoration_type="ENP", new_resid=nothing, new_segid=nothing, new_pdbname=nothing
    )
  
    decoration = isnothing(decoration) ? getPDB(decoration_type) : decoration
    tmpPDB = PDBTools.readPDB(decoration) ## loading the decoration PDB to make apply the transformations

    idxs, ref0_coord, ref1_coord, ini0_coord, ini1_coord = get_position_info(reference, decoration, decoration_type, resid=resid, segid=segid)

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
        tmpPDB[at].resnum = ifelse(isnothing(new_resid), resid, new_resid)
        tmpPDB[at].segname = ifelse(isnothing(new_segid), segid, new_segid)
    end

    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname
    PDBTools.writePDB(tmpPDB, new_pdbname)

    return new_pdbname

end






