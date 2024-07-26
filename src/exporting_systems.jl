"""

    _exporting_PDBfile(n::Int64, tmpfile_list::Vector{String}; topology_file="toppar/top_all36_carb.rtf", covalent=true, vmd="vmd")

This

## Arguments

- `n::Int64`: is the number of β-Glc units in each cellulose chain.
- `tmpfile_list::Vector{String}`: is the vector with all PDB files for each cellulose chain.

### Examples

```jldoctest

julia > transformingPBC(59, 5, 7)

```

"""

function _exporting_PDBfile(n::Int64, tmpfile_list::Vector{String}; phase="Iβ", covalent=true, vmd="vmd", topology_file="toppar/top_all36_carb.rtf")

    fraglist = collect(1:1:length(tmpfile_list))
    monomer, prev_monomer = 0, 0;

    pdb = joinpath(tempdir(), "cellulose.pdb")
    psf = joinpath(tempdir(), "cellulose.psf")
    tcl = joinpath(tempdir(), "cellulose.tcl")

            
    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput, "package require psfgen \n")
    Base.write(vmdinput, "topology $(topology_file) \n")
    Base.write(vmdinput, " \n")
    for frag in fraglist
        segname = "M$frag"
        Base.write(vmdinput, "segment $segname { \n")
        Base.write(vmdinput, "    pdb $(tmpfile_list[frag]) \n")
        Base.write(vmdinput, "} \n")
        monomer = n
        while monomer > 1
            prev_monomer = monomer - 1
            Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$prev_monomer \n")
            monomer -= 1
            if monomer == 1 && covalent == true
                Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$n \n")
            end
        end
    end
    Base.write(vmdinput, " \n")
    Base.write(vmdinput, "regenerate angles dihedrals \n")
    Base.write(vmdinput, " \n")
    for frag in fraglist
        segname = "M$frag"
        Base.write(vmdinput, "coordpdb $(tmpfile_list[frag]) $segname \n")
    end
    Base.write(vmdinput, " \n")
    Base.write(vmdinput, "guesscoord \n")
    Base.write(vmdinput, " \n")
    Base.write(vmdinput, "writepsf $(psf) \n")
    Base.write(vmdinput, "writepdb $(pdb) \n")
    Base.close(vmdinput)
    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")


    return vmdoutput
end