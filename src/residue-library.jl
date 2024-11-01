##struct Patch
##    symbol::Symbol
##    name::String
##    ref_atom_0::String
##    ref_atom_1::String
##    new_atom_0::String
##    new_atom_1::String
##end
##
##Base.Broadcast.broadcastable(patch::Patch) = Ref(patch)
##const patchs = Dict{String, Patch}([
##    ["ENP", "EthNP", "pEtN", "phosphoethanolamine"] .=> Patch(:PCEL, "PCEL", "O6", "HO6", "P", "O1");
##])

"""
    get_position_info(reference::String, decoration::String, residue::String; resid=nothing, segid=nothing)

    Takes the `main chain` (pdb1) and `decoration` (pdb2) PDB files and returns the positional data needed to translate and to rotate the `residue` before linkage.
    You must include the PDB filename of the decoration as argument for this method.

    # Arguments
    - `reference::String`: The name of the PDB file of the chain.
    - `decoration::String`: The name of the PDB file of the decoration.
    - `residue::String`: The name of the decoration residue. Currently, only `ENP`.
    - `resid::Int64`: The residue number of the main chain monomer that you will be modified. The default is: if `nothing`, the first residue number.
    - `segid::String`: The segment name of the main chain monomer that you will be modified. The default is: if `nothing`, the first segment name.
"""
function get_position_info(reference::String, decoration::String, residue::String; resid=nothing, segid=nothing)

    pdb1 = PDBTools.readPDB(reference)
    pdb2 = PDBTools.readPDB(decoration)

    resid = isnothing(resid) ? pdb1[1].resnum : resid
    segid = isnothing(segid) ? pdb1[1].segname : segid

    if residue == "ENP" || residue == "phosphoethanolamine" || residue == "EthNP" || residue == "pEtN"
        idx1 = PDBTools.index.(
                PDBTools.select(pdb1, (atom -> atom.name == "O6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx2 = PDBTools.index.(
                PDBTools.select(pdb1, (atom -> atom.name == "HO6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx3 = PDBTools.index.(
                PDBTools.select(pdb2, (atom -> atom.name == "P"))
            )[1]
        
        idx4 = PDBTools.index.(
                PDBTools.select(pdb2, (atom -> atom.name == "O1P"))
            )[1]
    end

    ## main
    coord1 = PDBTools.coor.(pdb1)[idx1]
    coord2 = PDBTools.coor.(pdb1)[idx2]
    ## decoration
    coord3 = PDBTools.coor.(pdb2)[idx3]
    coord4 = PDBTools.coor.(pdb2)[idx4]

    return [idx1, idx2, idx3, idx4], coord1, coord2, coord3, coord4

end

"""
    get_position_info(reference::String, residue::String; resid=nothing, segid=nothing)

    Takes the `main chain` (pdb1) and `decoration` (pdb2) PDB files and returns the positional data needed to translate and to rotate the `residue` before linkage.
    You must not include the PDB filename of the decoration as argument for this method, the function will pick using `getPDB()`.

    # Arguments
    - `reference::String`: The name of the PDB file of the chain.
    - `residue::String`: The name of the decoration residue. Currently, only `ENP`.
    - `resid::Int64`: The residue number of the main chain monomer that you will be modified. The default is: if `nothing`, the first residue number.
    - `segid::String`: The segment name of the main chain monomer that you will be modified. The default is: if `nothing`, the first segment name.
"""
function get_position_info(reference::String, residue::String; resid=nothing, segid=nothing)

    pdb1 = PDBTools.readPDB(reference)
    pdb2 = PDBTools.readPDB(getPDB(residue))

    resid = isnothing(resid) ? pdb1[1].resnum : resid
    segid = isnothing(segid) ? pdb1[1].segname : segid

    if residue == "ENP" || residue == "phosphoethanolamine" || residue == "EthNP" || residue == "pEtN"
        idx1 = PDBTools.index.(
                PDBTools.select(pdb1, (atom -> atom.name == "O6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx2 = PDBTools.index.(
                PDBTools.select(pdb1, (atom -> atom.name == "HO6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx3 = PDBTools.index.(
                PDBTools.select(pdb2, (atom -> atom.name == "P"))
            )[1]
        
        idx4 = PDBTools.index.(
                PDBTools.select(pdb2, (atom -> atom.name == "O1P"))
            )[1]
    end

    ## main
    coord1 = PDBTools.coor.(pdb1)[idx1]
    coord2 = PDBTools.coor.(pdb1)[idx2]
    ## decoration
    coord3 = PDBTools.coor.(pdb2)[idx3]
    coord4 = PDBTools.coor.(pdb2)[idx4]

    return [idx1, idx2, idx3, idx4], coord1, coord2, coord3, coord4

end


"""
    getPDB(residue::String; pdbname=nothing)
    
    This function creates a PDB file for a given residue.

    # Arguments
    - `residue`: The residue for which the PDB file is to be created. Currently, only `ENP `.
    - `pdbname`: The name of the PDB file to be created. If `nothing`, it will generate a temporary file
"""
function getPDB(residue::String; pdbname=nothing)
    
    if isnothing(pdbname)
        pdbname = tempname() * ".pdb"
    end

    pdb = Base.open(pdbname, "w")

    Base.write(pdb, "HEADER    GENERATED BY getPDB()             2024-10-11  XXXX\n")
    if residue == "ENP" || residue == "phosphoethanolamine" || residue == "EthNP" || residue == "pEtN"
        Base.write(pdb, "HETATM    1  N   ENP     1       2.249  -1.408  -0.019  1.00  0.00      U    N\n")
        Base.write(pdb, "HETATM    2  HN1 ENP     1       2.708  -2.325  -0.058  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM    3  HN2 ENP     1       2.550  -0.931   0.841  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM    4  HN3 ENP     1       2.584  -0.845  -0.811  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM    5  C12 ENP     1       0.731  -1.493  -0.061  1.00  0.00      U    C\n")
        Base.write(pdb, "HETATM    6  H2A ENP     1       0.458  -1.998  -0.988  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM    7  H2B ENP     1       0.410  -2.084   0.797  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM    8  C11 ENP     1       0.169  -0.086  -0.022  1.00  0.00      U    C\n")
        Base.write(pdb, "HETATM    9  H1A ENP     1       0.576   0.524  -0.843  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM   10  H1B ENP     1       0.441   0.404   0.928  1.00  0.00      U    H\n")
        Base.write(pdb, "HETATM   11  P   ENP     1      -2.051   1.236  -0.086  1.00  0.00      U    P\n")
        Base.write(pdb, "HETATM   12  O4P ENP     1      -1.143   2.237  -0.755  1.00  0.00      U    O\n")
        Base.write(pdb, "HETATM   13  O3P ENP     1      -3.455   0.906  -0.453  1.00  0.00      U    O\n")
        Base.write(pdb, "HETATM   14  O2P ENP     1      -1.207  -0.267  -0.142  1.00  0.00      U    O\n")
        Base.write(pdb, "HETATM   15  O1P ENP     1      -1.865   1.459   1.564  0.00  0.00      U    O\n")
        Base.write(pdb, "HETATM   16  HOP ENP     1      -1.632   2.399   1.669  0.00  0.00      U    H\n")
    else
        throw(ArgumentError("The residue $residue is still not covered by this script."))
    end
    Base.write(pdb, "END\n")

    Base.close(pdb)

    return pdbname

end