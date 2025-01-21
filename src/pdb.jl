#1         2         3         4         5         6         7         8
#12345678901234567890123456789012345678901234567890123456789012345678901234567890
#HETATM 1357 MG    MG   168       4.669  34.118  19.123  1.00  3.16          MG2+
#HETATM 3835 FE   HEM     1      17.140   3.115  15.066  1.00 14.14          FE3+
#HETATM 1407  CA  BLE P   1      14.625  32.240  14.151  1.09 16.76           C

#ATOM      1  C       X   1       0.110   0.434   3.923  0.00  0.00           C
#HETATM    1  C1  BGC     1       2.515   8.013   3.923  1.00  0.00           C         after the treatment

"""
    cleanPDB!(pdbname::String, atomnames::Vector{String})

Clean the PDB file by changing the atomnames, resnames, resid, and chain. The crystallographic coordinates are marked as 1.00.
"""
function cleanPDB!(pdbname::String, atomnames::Vector{String})

    pdbfile = String[]

    iatom, resid = 1, 1

    push!(pdbfile, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1")
    for line in split(read(pdbname, String), "\n")

        if !startswith(line, "ATOM  ")
            continue
        end

        newline = "HETATM" * line[7:end]
        newline = string(newline[1:13], rpad(atomnames[iatom], 3), newline[17:end])         
        newline = string(newline[1:17], "BGC", newline[21:end])                             
        newline = string(newline[1:21], " ", newline[23:end])                               
        newline = string(newline[1:22], lpad(resid, 4), newline[27:end])                    
        newline = string(newline[1:56], "1.00", newline[61:end])                            

        push!(pdbfile, newline)

        iatom += 1
        if iatom > length(atomnames)
            iatom = 1
            resid += 1
        end
    end
    push!(pdbfile, "END")

    open(pdbname, "w") do file
        for line in pdbfile
            println(file, line)
        end
    end

    return nothing
end


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

    vmdoutput = execVMD(vmd, tcl)
    
    return vmdDebug ? vmdoutput : pdbnames
end

"""
    writePDB(monomers::Int64, pdbnames::Vector{String}; phase="Iβ", covalent=true, check_inversion=false, topology_file=DEFAULT_CARB_TOPOLOGY_FILE, vmd="vmd", vmdDebug=true)

Writes the final PSF and PDB files from the cellulose fragments.
"""
function writePDB(
                monomers::Int64,
                pdbnames::Vector{String};
                pdbname=nothing,
                phase="Iβ", covalent=true, hasInversion=false,
                topology_file=DEFAULT_CARB_TOPOLOGY_FILE,
                vmd="vmd", vmdDebug=true
            )
    
    pdb = isnothing(pdbname) ? tempname() * ".pdb" : pdbname
    psf, tcl = replace(pdb, ".pdb" => ".psf"), replace(pdb, ".pdb" => ".tcl")

    itoken = hasInversion ? isinverted.(pdbnames) : nothing
            
    Base.open(tcl, "w") do file
        println(file, """
        package require psfgen
        topology $topology_file
        """)

        for (id, pdbname) in enumerate(pdbnames)
            
            println(file, """
            segment M$id { 
                pdb $pdbname 
            } 
            """)

            inverting = !isnothing(itoken) ? itoken[id] : false
            printching!(file, id=id, nresids=monomers, invert=inverting, phase=phase, covalent=covalent)

            println(file, """
            """) ## to be more clean on TCL
        end

        println(file, """
        regenerate angles dihedrals
        """)

        for (id, pdbname) in enumerate(pdbnames)
            println(file, "coordpdb $pdbname M$id")
        end
        
        println(file, """
        guesscoord

        writepsf $psf
        writepdb $pdb
        exit
        """)
    
    end
    
    vmdoutput = execVMD(vmd, tcl)
    
    return vmdDebug ? vmdoutput : (pdb, psf)
end