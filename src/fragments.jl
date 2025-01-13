"""
    picking_fragments(vmdoutput)

Picking the number of fragments (e.g. cellulose chains) inside the structure loaded on VMD. Usualy,
the VMD output has the type `Vector{SubString{String}}` such as:

### Examples

```jldoctest

julia > picking_fragments(vmdoutput)
julia > 59

```

"""
function picking_fragments(vmdoutput::Vector{SubString{String}})

    nfragments = nothing
    
    for id in eachindex(vmdoutput)
        str = string(vmdoutput[id])
        if occursin("Fragments:", str)
            nfragments = parse(Int64, split(str)[3])
        elseif isnothing(nfragments) && id == length(vmdoutput)
            error("The VMD output file does not contain the number of fragments.")
        end
    end
    
    return nfragments
    
end

"""

    writeXYZ(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    writeXYZ(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    writeXYZ(atoms::Vector{String}, xyzcoords::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)

This function is able to save the atomic information (labels + cartesian coordinates) on XYZ file.

## Arguments

- `atoms (::Vector{Vector{String}}, ::Vector{String})`: atoms names of the systems
- `x, y, z (::Vector{Vector{Float64}}, ::Vector{Float64})`: cartesian coordinates of the system
- `xyzcoords::Vector{Vector{Float64}}`: (x,y,z) of the cartesian coordinates

### Examples

```jldoctest

julia > writeXYZ(atoms, x, y, z)
julia > writeXYZ(atoms, xyz)

```

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

function rawXYZ(xyz::XYZs, dim::Vector{Int64}; phase="Iβ", exporting=true, natoms=nothing, xyzfile=nothing, vmd="vmd", output=false)
    
    xyz_trimmed = _trimming_xy(xyz, dim, phase=phase)
    xyz_expanded = _expanding_z(xyz_trimmed, dim[3], phase=phase)

    if !exporting
        return xyz_expanded
    else
        return writeXYZ(
                    xyz_expanded.atoms, xyz_expanded.x, xyz_expanded.y, xyz_expanded.z;
                    natoms=natoms, xyzfile=xyzfile, vmd=vmd, output=output
                )
    end
end