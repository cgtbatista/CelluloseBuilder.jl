"""

    _exporting_XYZfile(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    _exporting_XYZfile(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    _exporting_XYZfile(atoms::Vector{String}, xyzcoords::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)

This function is able to save the atomic information (labels + cartesian coordinates) on XYZ file.

## Arguments

- `atoms (::Vector{Vector{String}}, ::Vector{String})`: atoms names of the systems
- `x, y, z (::Vector{Vector{Float64}}, ::Vector{Float64})`: cartesian coordinates of the system
- `xyzcoords::Vector{Vector{Float64}}`: (x,y,z) of the cartesian coordinates

### Examples

```jldoctest

julia > _exporting_XYZfile(atoms, x, y, z)
julia > _exporting_XYZfile(atoms, xyz)

```

"""

function _exporting_XYZfile(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    newatoms = reduce(vcat, atoms)
    newx = reduce(vcat, x); newy = reduce(vcat, y); newz = reduce(vcat, z);
    return _exporting_XYZfile(newatoms, newx, newy, newz, vmd=vmd, natoms=natoms, xyzfile=xyzfile)
end

function _exporting_XYZfile(atoms::Vector{String}, xyzcoords::Vector{Vector{Float64}}; vmd="vmd", natoms=nothing, xyzfile=nothing)
    x, y, z = Float64[], Float64[], Float64[]
    for xyz in xyzcoords
        push!(x, xyz[1]); push!(y, xyz[2]); push!(z, xyz[3])
    end
    return _exporting_XYZfile(atoms, x, y, z, vmd=vmd, natoms=natoms, xyzfile=xyzfile)
end

function _exporting_XYZfile(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}; vmd="vmd", natoms=nothing, xyzfile=nothing)

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
     
    return xyzfile, vmdoutput

end



"""

    _PDBfragment_cleaning(atoms::Vector{String}, pdbfile::String, new_pdbfile::String)

Clean the raw PDB file and export a new one with the right configuration needed for CHARMM force field.

## Arguments

- `atoms::Vector{String}`: atoms names of the systems
- `x, y, z (::Vector{Vector{Float64}}, ::Vector{Float64})`: cartesian coordinates of the system
- `xyzcoords::Vector{Vector{Float64}}`: (x,y,z) of the cartesian coordinates

### Examples

```jldoctest

julia > _PDBfragment_cleaning(atoms, "/tmp/old.pdb", "/tmp/new.pdb")

```

"""

function _PDBfragment_cleaning(atoms::Vector{String}, pdbfile::String, new_pdbfile::String)

    ith_atom = 1; resid = 1;

    pdbdata = Base.split(Base.read(pdbfile, String), "\n")
    pdb = open(new_pdbfile, "w")

    Base.write(pdb, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")

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
        new_irow = replace(new_irow, _pdbcolumn_resname => "BGLC  ")
        # editing residue sequence numbers
        _pdbcolumn_resnum = SubString(new_irow, 26-length("$resid")+1 : 26)
        new_irow = replace(new_irow, _pdbcolumn_resnum*"    " => "$resid"*"    ") ## This dummy space is escential to not replace all numbers (x coords begin on 31 string position)
        # editing the segname
        _pdbcolumn_segname = SubString(new_irow, findfirst("  0.00  0.00           ", irow))
        new_irow = replace(new_irow, _pdbcolumn_segname => "  1.00  0.00           ")
        
        # writting the new configuration
        Base.write(pdb, "$new_irow\n")
        
        ith_atom += 1
        if ith_atom > length(atoms)
            ith_atom = 1
            resid += 1
        end
    end
    Base.write(pdb, "END")
    
    Base.close(pdb)

    return nothing

end


"""

    _XYZfragments_2_PDB(xyzfile::String, selection::String; vmd="vmd")

Convert each fragment on XYZ file on an unique raw PDB file. The filename will be exported, so there is
a way to edit the PDB data.

## Arguments

- `xyzfile::String`: the XYZ filename.
- `vmdselection::String`: the VMD selection type.

### Examples

```jldoctest

julia > _XYZfragments_2_PDB("system.xyz", fragment_list)

```

"""

function _XYZfragments_2_PDB(xyzfile::String, selection::String; vmd="vmd")
    
    filename = tempname()

    vmdinput_file = tempname() * ".tcl"
    vmdinput = open(vmdinput_file, "w")

    Base.write(vmdinput, "mol new \"$xyzfile\" \n")
    units = split(selection, " ")
    for u in units
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $u\" ] \n")
        Base.write(vmdinput, "\$sel writepdb \"$(filename * "_" * u * ".pdb")\" \n")
    end
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)
    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    return vmdoutput, filename

end