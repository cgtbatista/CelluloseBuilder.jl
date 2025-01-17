function cleanPDB!(pdbname::String, atoms::Vector{String})

    pdbdata = Base.split(Base.read(pdbname, String), "\n")

    new_pdbdata = Vector{String}()

    push!(new_pdbdata, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1")

    ith_atom, resid = 1, 1
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
        new_irow = replace(new_irow, _pdbcolumn_resname => "BGC   ")
        # editing residue sequence numbers
        _pdbcolumn_resnum = SubString(new_irow, 26-length("$resid")+1 : 26)
        new_irow = replace(new_irow, _pdbcolumn_resnum*"    " => "$resid"*"    ")
        # editing the segname
        _pdbcolumn_segname = SubString(new_irow, findfirst("  0.00  0.00           ", irow))
        new_irow = replace(new_irow, _pdbcolumn_segname => "  1.00  0.00           ")
        
        # writting the new configuration
        push!(new_pdbdata, "$new_irow")
        
        ith_atom += 1
        if ith_atom > length(atoms)
            ith_atom = 1
            resid += 1
        end
    end
    push!(new_pdbdata, "END")

    open(pdbname, "w") do pdb
        for irow in new_pdbdata
            println(pdb, irow)
        end
    end    

    return nothing
end

#HETATM    1  C1  BGC     1      11.074  40.399   0.449  1.00  0.00           C
#HETATM    1 C1   BGC     1  1      11.074  40.399   0.  1.00  0.00               C


# function cleanPDB!(pdbname::String, atoms::Vector{String})

#     pdbdata = split(read(pdbname, String), "\n")

#     new_pdbdata = String[]
#     push!(new_pdbdata, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1")

#     ith_atom, resid = 1, 1

#     for irow in pdbdata

#         if !startswith(irow, "ATOM  ")
#             continue
#         end

#         new_irow = "HETATM" * irow[7:end]                                               # updating atom pattern to HETATM
#         atomname = rpad(atoms[ith_atom], 4)                                             # granting right padding
#         new_irow = string(new_irow[1:12], atomname, new_irow[17:end])
#         new_irow = string(new_irow[1:17], "BGC   ", new_irow[21:end])                   # updating residue name
#         new_irow = string(new_irow[1:22], lpad(resid, 4), new_irow[27:end])             # updating residue number
#         new_irow = string(new_irow[1:54], "  1.00  0.00           ", new_irow[77:end])

#         push!(new_pdbdata, new_irow)

#         ith_atom += 1
#         if ith_atom > length(atoms)
#             ith_atom = 1
#             resid += 1
#         end
#     end

#     push!(new_pdbdata, "END")

#     open(pdbname, "w") do pdb
#         for line in new_pdbdata
#             println(pdb, line)
#         end
#     end

#     return nothing
# end


"""
    fragPDBs(xyzfile::String, selection::String; vmd="vmd")

Convert each fragment on XYZ file on an unique raw PDB file. The filename will be exported, so there is a way to edit the PDB data.
"""
function fragPDBs(xyzname::String; fragments=nothing, vmd="vmd", vmdDebug=false)

    basename = tempname()

    fragments = isnothing(fragments) ? Base.OneTo(picking_fragments(xyzname, vmd=vmd)) : fragments
    pdbnames = [ "$(basename)_$(frag).pdb" for frag in fragments ]

    tcl = basename * ".tcl"
    open(tcl, "w") do file
        println(file, "mol new \"$xyzname\"")
        for (f, pdbname) in enumerate(pdbnames)
            println(file, """
            set sel [atomselect top "fragment $(f-1)"]
            \$sel writepdb "$pdbname"
            """)
        end
        println(file, "exit")
    end

    vmdoutput = split(read(`$vmd -dispdev text -e $tcl`, String), "\n")
    
    return vmdDebug ? vmdoutput : pdbnames
end
