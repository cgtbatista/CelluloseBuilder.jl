"""
    writePDB(monomers::Int64, pdbnames::Vector{String}; phase="Iβ", covalent=true, check_inversion=false, topology_file=DEFAULT_CARB_TOPOLOGY_FILE, vmd="vmd", vmdDebug=true)

Writes the final PSF and PDB files from the cellulose fragments.
"""
function writePDB(
                monomers::Int64,
                pdbnames::Vector{String};
                phase="Iβ", covalent=true, check_inversion=false,
                topology_file=DEFAULT_CARB_TOPOLOGY_FILE,
                vmd="vmd", vmdDebug=true
            )

    pdb, psf, tcl = joinpath(tempdir(), "cellulose.pdb"), joinpath(tempdir(), "cellulose.psf"), joinpath(tempdir(), "cellulose.tcl")

    itoken = check_inversion ? isinverted.(pdbnames) : nothing
            
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
            """) ## to be more clean on the tcl file
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
    
    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")
    if vmdDebug
        return vmdoutput
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