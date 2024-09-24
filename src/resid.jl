"""
    get_residuePDB(residue::String; pdbname=nothing)
    
    This function creates a PDB file for a given residue.

#### Arguments
- `residue`: The residue for which the PDB file is to be created. Currently, only `PETN`.
- `pdbname`: The name of the PDB file to be created. If `nothing`, it will generate a temporary file
"""
function get_residuePDB(residue::String; pdbname=nothing)
    
    if isnothing(pdbname)
        pdbname = tempname() * ".pdb"
    end

    pdb = Base.open(pdbname, "w")

    Base.write(pdb, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")

    if residue == "PETN"
        Base.write(pdb, "ATOM      1  N   PETNU   1       6.758  -1.101   4.264  1.00  0.00      U    N\n")
        Base.write(pdb, "ATOM      2  HN1 PETNU   1       7.717  -0.971   3.936  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      3  HN2 PETNU   1       6.637  -0.657   5.194  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      4  HN3 PETNU   1       6.023  -0.537   3.647  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      5  C2  PETNU   1       6.341  -2.548   4.259  1.00  0.00      U    C\n")
        Base.write(pdb, "ATOM      6  H2A PETNU   1       7.200  -3.148   3.946  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      7  H2B PETNU   1       6.059  -2.824   5.275  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      8  C1  PETNU   1       5.171  -2.812   3.308  1.00  0.00      U    C\n")
        Base.write(pdb, "ATOM      9  H1A PETNU   1       5.096  -3.895   3.175  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM     10  H1B PETNU   1       5.368  -2.353   2.333  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM     11  P   PETNU   1       3.704  -0.781   4.106  1.00  0.00      U    P\n")
        Base.write(pdb, "ATOM     12  O3  PETNU   1       2.278  -0.388   4.235  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     13  O4  PETNU   1       4.698  -0.069   3.168  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     14  O1  PETNU   1       3.959  -0.457   5.017  0.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     15  O2  PETNU   1       3.901  -2.400   3.806  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     16  HO  PETNU   1       4.785  -0.937   5.314  0.00  0.00      U    H\n")
    else
        throw(ArgumentError("The residue $residue is still not covered by this script."))
    end

    Base.write(pdb, "END\n")

    Base.close(pdb)

    return pdbname

end

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
    ith_patch = maximum(PDBTools.resnum.(pdb))
    for i in resid
        ith_patch += ith_patch
        new_patch = matching_residue(chain_pdbname, i, chain_name, segid, decoration=decoration, new_resid=ith_patch, new_chain_name=new_chain_name, new_segid=new_segid)
        push!(patchings, [i, ith_patch])
        push!(pdb_decorations, new_patch)
    end

    patches_pdb = combinePDBs(pdb_decorations, new_pdbfile=patches_pdbname)

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

function combinePDBs(pdbfiles::Vector{String}; new_pdbfile=nothing)

    new_pdbfile = isnothing(new_pdbfile) ? tempname() * ".pdb" : new_pdbfile
    pdb = Base.open(new_pdbfile, "w")

    Base.write(pdb, "REMARK  1 PDBCOMBINE\n")

    atom_index = 1

    for pdbfile in pdbfiles

        tmpPDB = PDBTools.readPDB(pdbfile)
        for i in eachindex(tmpPDB)
            tmpPDB[i].index = atom_index
            tmpPDB[i].index_pdb = atom_index
            atom_index += 1
        end

        PDBTools.writePDB(tmpPDB, pdbfile)

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

function matching_residue(chain_pdbname::String, resid::Int64, chain_name::String, segid::String; decoration="PETN",
        new_resid=nothing, new_chain_name=nothing, new_segid=nothing,
        pdbname=nothing, new_pdbname=nothing
    )
  
    res_pdbname = isnothing(pdbname) ? get_residuePDB(decoration) : pdbname
    res = PDBTools.readPDB(res_pdbname)  

    idxs, ref0_coord, ref1_coord, ini0_coord, ini1_coord = decoration_library(chain_pdbname, res_pdbname, decoration, resid=resid, segid=segid)

    # Rotating the residue coordinates following the β-Glc residue orientation...
    v_residue = ini1_coord - ini0_coord
    v_target = ref0_coord - ref1_coord
    rotated_coords = rotate_residue(
            PDBTools.coor.(res),
            ref1_coord,
            R_matrix(v_target, v_residue)
        )
    # Translating the residue coordinates to the specify position...
    v_trans = ref0_coord - rotated_coords[idxs[3]]
    translated_coords = translate_residue(
            rotated_coords, v_trans
        )

    for at in eachindex(translated_coords)
        res[at].x = translated_coords[at][1]
        res[at].y = translated_coords[at][2]
        res[at].z = translated_coords[at][3]
        res[at].chain = ifelse(isnothing(new_chain_name), chain_name, new_chain_name)
        res[at].resnum = ifelse(isnothing(new_resid), resid, new_resid)
    end

    new_pdbname = isnothing(new_pdbname) ? tempname() * ".pdb" : new_pdbname
    PDBTools.writePDB(res, new_pdbname)

    vmdoutput = updating_segid(
            new_pdbname,
            ifelse(isnothing(new_segid), segid, new_segid),
            new_pdbname=new_pdbname
        ) ## não sei o porquê, mas não consigo salvar o nome do segmento no PDBTools: `res[at].segname = ifelse(isnothing(new_segid), segid, new_segid)`

    return new_pdbname

end

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

function honeycomb_assembly(center::Vector{Float64}, r::Float64; space=12.)

    centers = Float64[]

    # Coordenadas da fibrila central
    x0, y0, z0 = center[1], center[2], center[3]

    # Lista para armazenar os centros das 7 fibrilas
    centers = []

    adj_space = space + 2 * r
    
    push!(centers, [x0, y0, z0])
    push!(centers, [x0 - adj_space, y0, z0])
    push!(centers, [x0 + adj_space, y0, z0])
    push!(centers, [x0 - adj_space/2, y0 + adj_space * √3/2, z0])
    push!(centers, [x0 + adj_space/2, y0 + adj_space * √3/2, z0])
    push!(centers, [x0 - adj_space/2, y0 - adj_space * √3/2, z0])
    push!(centers, [x0 + adj_space/2, y0 - adj_space * √3/2, z0])

    return centers

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

"""
    decoration_library(chain_pdbname, decoration_pdbname, decoration::String)

    This function is a library of decorations that can be added to the cellulose chain. The decorations are: PETN...
"""
function decoration_library(chain_pdbname, decoration_pdbname, decoration::String; resid=1, segid=nothing)

    main_chain = PDBTools.readPDB(chain_pdbname)
    side_resid = PDBTools.readPDB(decoration_pdbname)

    segid = isnothing(segid) ? main_chain[1].segname : segid

    if decoration == "PETN"
        idx1 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "O6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx2 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "HO6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx3 = PDBTools.index.(
                PDBTools.select(side_resid, by = (atom -> atom.name == "P"))
            )[1]
        
        idx4 = PDBTools.index.(
                PDBTools.select(side_resid, by = (atom -> atom.name == "O1"))
            )[1]
    end

    ## main
    coord1 = PDBTools.coor.(main_chain)[idx1]
    coord2 = PDBTools.coor.(main_chain)[idx2]
    ## decoration
    coord3 = PDBTools.coor.(side_resid)[idx3]
    coord4 = PDBTools.coor.(side_resid)[idx4]

    return [idx1, idx2, idx3, idx4], coord1, coord2, coord3, coord4

end

"""
    R_matrix(v1, v2)

    Picks two vectors and gives the rotation matrix `R` needed to rotate the vector v2 on v1 (the mirror vector).
    If the ```||v_cross|| ≈ 0```, we do not need to rotate the residue. Otherwise, we will pick a generalize rotation matrix `R`.
"""
function R_matrix(v1, v2)

    # Normalizing the vectors
    v1 /= norm(v1)
    v2 /= norm(v2)

    # cross rotation vector
    crossvector = cross(v1, v2)
    θ = acos(
            clamp(dot(v1, v2), -1.0, 1.0)
        ) # just to ensure the angle domain


    # checking with there is need in computation of the cross matrix
    if norm(crossvector) ≈ 0
        return I(3) # return the identity
    else
        crossvector /= norm(crossvector)
    end

    # cross matrix
    CrossMatrix = [ 0 -crossvector[3] crossvector[2]; crossvector[3] 0 -crossvector[1]; -crossvector[2] crossvector[1] 0 ]
    R = I(3) + sin(θ) * CrossMatrix + (1-cos(θ)) * (CrossMatrix * CrossMatrix)

    return R

end

"""
    translate_residue(coords, translation_vector)

    Translate the residue coordinates by a given directional vector.
"""
function translate_residue(coords, translation_vector)

    for xyz in eachindex(coords)
        coords[xyz] += translation_vector
    end

    return coords
end

"""
    rotate_residue(coords, reference, rotation_matrix)

    Rotate the residue coordinates by a given rotation matrix. The `ref_center` is a reference point that we will lock while rotating the other atoms.
"""
function rotate_residue(coords, ref_center, rotation_matrix)

    for xyz in eachindex(coords)
        coords[xyz] = ref_center + rotation_matrix * (coords[xyz] - ref_center)
    end

    return coords
end

