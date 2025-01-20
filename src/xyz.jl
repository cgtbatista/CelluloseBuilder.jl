function rawXYZ(crystal::CrystalXYZ, xyzsizes::Vector{Int64}; phase="Iβ", exporting=true, natoms=nothing, xyzfile=nothing, vmd="vmd", output=false)
    
    atemp, xtemp, ytemp, ztemp = xy_pruning(crystal.atoms, crystal.x, crystal.y, crystal.z, xyzsizes, phase)
    atoms, x, y, z = z_expansion(atemp, xtemp, ytemp, ztemp, xyzsizes, phase)

    if !exporting
        return atoms, x, y, z
    else
        return writeXYZ(
                    atoms, x, y, z;
                    natoms=natoms, xyzfile=xyzfile, vmd=vmd, output=output
                )
    end
end

function pbcXYZ(
            xyzname::String,
            xyzsizes::Vector{Int64};
            phase="Iβ", structure=nothing, fibril=nothing, new_xyzname=nothing,
            vmd="vmd", vmdDebug=false
        )
    
    new_xyzname = isnothing(new_xyzname) ? tempname() * ".xyz" : new_xyzname
    
    structure = isnothing(structure) ? "fibril" : structure
    
    valid_structure = in(lowercase(structure), Set(["fibril", "f", "monolayer", "m", "center", "c", "origin", "o"]))
    if !valid_structure
        error("The structure $structure is not valid.")
    end
    
    isfibril = in(lowercase(structure), Set(["fibril", "f"]))

    if isfibril
        fibril = isnothing(fibril) ? "34566543" : fibril

        key = (lowercase(phase), lowercase(fibril))
        rchains = if haskey(microfibril, key)
                parse.(Int64, split(microfibril[key], " "))
            else
                error("The fibril $fibril is not implemented yet on the microfibril dictionary.")
        end
    end
   
    if !isfibril
        rchains = monolayer(xyzsizes, phase=phase, structure=structure)
    end

    tcl = tempname() * ".tcl"

    vmdinput = open(tcl, "w")
    Base.write(vmdinput, """
    mol new $xyzname
    set sel [ atomselect top \"fragment $(join(rchains, " "))\" ]
    \$sel writexyz $new_xyzname
    exit
    """)
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if vmdDebug
        return vmdoutput
    else
        return new_xyzname, rchains
    end
end

"""
    pbcXYZ(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")

Applies the periodic boundary conditions on the cellulose crystal over *xy-plane* using the number of fragments (e.g. chains) as parameter.
It returns cellulose **XYZ file** needed to prepare the PDB and PSF files.
"""
function pbcXYZ(
            xyzname::String,
            fragments::Int64,
            xyzsizes::Vector{Int64};
            phase="Iβ", pbc=nothing, new_xyzname=nothing,
            vmd="vmd", vmdDebug=false
        )    
    
    valid_phase = in(lowercase(phase), Set(["ib", "iβ", "ii", "ia", "iα", "iii", "iii_i", "iii_i", "iii_i"]))
    if !valid_phase
        error("The phase $phase is not implemented yet.")
    end

    new_xyzname = isnothing(new_xyzname) ? tempname() * ".xyz" : new_xyzname
        
    if in(lowercase(phase), Set(["ib", "iβ", "ii"]))        
        fchains = forbbiden(fragments, xyzsizes, pbc=pbc)
        rchains = setdiff(0:fragments-1, fchains)
    end

    if in(lowercase(phase), Set(["ia", "iα", "iii", "iii_i", "iiii"]))
        rchains = 0:fragments-1
    end

    tcl = tempname() * ".tcl"

    open(tcl, "w") do file
        println(file, """
        mol new "$xyzname"
        set sel [ atomselect top "fragment $(join(rchains, " "))" ]
        \$sel writexyz $new_xyzname
        exit
        """)
    end
    
    vmdoutput = execVMD(vmd, tcl)

    if vmdDebug
        return vmdoutput
    else
        return new_xyzname, rchains
    end
end


"""

    writeXYZ(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    writeXYZ(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    writeXYZ(atoms::Vector{String}, xyzcoords::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)

This function is able to save the atomic information (labels + cartesian coordinates) on XYZ file.
"""
function writeXYZ(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    newatoms = reduce(vcat, atoms)
    newx = reduce(vcat, x); newy = reduce(vcat, y); newz = reduce(vcat, z);
    return writeXYZ(newatoms, newx, newy, newz, vmd=vmd, natoms=natoms, xyzfile=xyzfile)
end

function writeXYZ(atoms::Vector{String}, xyzcoords::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    x, y, z = Float64[], Float64[], Float64[]
    for xyz in xyzcoords
        push!(x, xyz[1]); push!(y, xyz[2]); push!(z, xyz[3])
    end
    return writeXYZ(atoms, x, y, z, vmd=vmd, natoms=natoms, xyzfile=xyzfile)
end

function writeXYZ(
                atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
                natoms=nothing, xyzfile=nothing, vmd="vmd", output=false
            )

    if isnothing(natoms); natoms = length(atoms); end
    if isnothing(xyzfile); xyzfile = tempname() * ".xyz"; end
    
    structure = open(xyzfile, "w")
    Base.write(structure, "$(length(atoms)) 1\n")
    Base.write(structure, " \n")
    for (at, i, j, k) in zip(atoms, x, y, z)
        Base.write(structure, " $at $i $j $k\n")
    end
    Base.write(structure, " \n")
    Base.close(structure)

    vmdinput_file = tempname() * ".tcl"
    vmdinput = open(vmdinput_file, "w")
    Base.write(vmdinput, "mol new \"$xyzfile\" \n")
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    if output
        return vmdoutput
    else
        return xyzfile, picking_fragments(vmdoutput)
    end

end