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
##    ["pEtN", "PETN", "petn"] .=> Patch(:PCEL, "PCEL", "O6", "HO6", "P", "O1");
##])

"""
    decoration_library(chain_pdbname, decoration_pdbname, decoration::String)

    This function is a library of decorations that can be added to the cellulose chain. The decorations are: PETN...
"""
function decoration_library(main_pdbname::String, decoration_pdbname::String, label::String; resid=nothing, segid=nothing)

    main_chain = PDBTools.readPDB(main_pdbname)
    decoration = PDBTools.readPDB(decoration_pdbname)

    resid = isnothing(resid) ? main_chain[1].resnum : resid
    segid = isnothing(segid) ? main_chain[1].segname : segid

    if label == "PETN"
        idx1 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "O6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx2 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "HO6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx3 = PDBTools.index.(
                PDBTools.select(decoration, by = (atom -> atom.name == "P"))
            )[1]
        
        idx4 = PDBTools.index.(
                PDBTools.select(decoration, by = (atom -> atom.name == "O1"))
            )[1]
    end

    ## main
    coord1 = PDBTools.coor.(main_chain)[idx1]
    coord2 = PDBTools.coor.(main_chain)[idx2]
    ## decoration
    coord3 = PDBTools.coor.(decoration)[idx3]
    coord4 = PDBTools.coor.(decoration)[idx4]

    return [idx1, idx2, idx3, idx4], coord1, coord2, coord3, coord4

end


function decoration_library(main_pdbname::String, label::String; resid=nothing, segid=nothing)

    main_chain = PDBTools.readPDB(main_pdbname)
    decoration = PDBTools.readPDB(getPDB(label))

    resid = isnothing(resid) ? main_chain[1].resnum : resid
    segid = isnothing(segid) ? main_chain[1].segname : segid

    if label == "PETN"
        idx1 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "O6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx2 = PDBTools.index.(
                PDBTools.select(main_chain, by = (atom -> atom.name == "HO6" && atom.resnum == resid && atom.segname == segid))
            )[1]

        idx3 = PDBTools.index.(
                PDBTools.select(decoration, by = (atom -> atom.name == "P"))
            )[1]
        
        idx4 = PDBTools.index.(
                PDBTools.select(decoration, by = (atom -> atom.name == "O1"))
            )[1]
    end

    ## main
    coord1 = PDBTools.coor.(main_chain)[idx1]
    coord2 = PDBTools.coor.(main_chain)[idx2]
    ## decoration
    coord3 = PDBTools.coor.(decoration)[idx3]
    coord4 = PDBTools.coor.(decoration)[idx4]

    return [idx1, idx2, idx3, idx4], coord1, coord2, coord3, coord4

end