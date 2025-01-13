"""

    _cleaning_PDBfragment(atoms::Vector{String}, pdbfile::String, new_pdbfile::String)

Clean the raw PDB file and export a new one with the right configuration needed for CHARMM force field.

## Arguments

- `atoms::Vector{String}`: atoms names of the systems
- `x, y, z (::Vector{Vector{Float64}}, ::Vector{Float64})`: cartesian coordinates of the system
- `xyzcoords::Vector{Vector{Float64}}`: (x,y,z) of the cartesian coordinates

### Examples

```jldoctest

julia > _cleaning_PDBfragment(atoms, "/tmp/old.pdb", "/tmp/new.pdb")

```

"""

function _cleaning_PDBfragment(atoms::Vector{String}, pdbfile::String, new_pdbfile::String)

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
        new_irow = replace(new_irow, _pdbcolumn_resname => "BGC   ")
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

function _XYZfragments_2_PDB(xyzfile::String, selection::String, units::Int64; vmd="vmd")
    
    filename = tempname()

    vmdinput_file = tempname() * ".tcl"
    vmdinput = open(vmdinput_file, "w")

    Base.write(vmdinput, "mol new \"$xyzfile\" \n")
    fragments = split(selection, " ")
    for u in collect(1:1:units)
        frag = fragments[u]
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $(u-1)\" ] \n")
        Base.write(vmdinput, "\$sel writepdb \"$(filename * "_" * frag * ".pdb")\" \n")
    end
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)
    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    return vmdoutput, filename

end