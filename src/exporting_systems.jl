function _check_inversion(filename::String; vmd="vmd")

    tcl = tempname() * ".tcl"
    inversion_value = 0

    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput, "mol new $filename \n")
    Base.write(vmdinput, " \n")
    Base.write(vmdinput, "set sel [atomselect top \"resid 1 and name C5 O5\"] \n")
    Base.write(vmdinput, "set diff_value [expr [lindex \$sel 1] - [lindex \$sel 0] ] \n")
    Base.write(vmdinput, "if { \$diff_value < 0 } { puts \"Invertion: 1\" } else { puts \"Invertion: 0\" } \n")
    Base.write(vmdinput, " \n")
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    for id in eachindex(vmdoutput)
        str = string(vmdoutput[id])
        if occursin("Inversion:", str)
            inversion_value = parse(Int64, split(str)[2])
        elseif isnothing(inversion_value) && id == length(vmdoutput)
            error("The VMD output file does not contain the inversion code.")
        end
    end

    return inversion_value

end


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

function _exporting_PDBfile(n::Int64, tmpfile_list::Vector{String}; phase="Iβ", covalent=true, vmd="vmd", topology_file=DEFAULT_CARB_TOPOLOGY_FILE, check_inversion=false)

    fraglist = collect(1:1:length(tmpfile_list))
    monomer, prev_monomer = 0, 0;

    pdb = joinpath(tempdir(), "cellulose.pdb")
    psf = joinpath(tempdir(), "cellulose.psf")
    tcl = joinpath(tempdir(), "cellulose.tcl")

    segid_token = zeros(Int64, length(tmpfile_list))
    if check_inversion
        for frag in fraglist
            segid_token[frag] = _check_inversion(tmpfile_list[frag], vmd=vmd)
        end
    end
            
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
        _inversion_patch_token = segid_token[frag] != 0
        if phase == "II" && _inversion_patch_token
            while monomer > 1
                prev_monomer = monomer - 1
                Base.write(vmdinput, "patch 14bb $segname:$prev_monomer $segname:$monomer \n")
                monomer -= 1
                if monomer == 1 && covalent == true
                    Base.write(vmdinput, "patch 14bb $segname:$n $segname:$monomer \n")
                end
            end
        else
            while monomer > 1
                prev_monomer = monomer - 1
                Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$prev_monomer \n")
                monomer -= 1
                if monomer == 1 && covalent == true
                    Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$n \n")
                end
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