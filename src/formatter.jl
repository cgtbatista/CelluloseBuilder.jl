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


function _XYZfragments_2_PDB(xyzfile::String, selection::String; vmd="vmd")
    
    rawname = tempname()

    vmdinput_file = tempname() * ".tcl"
    vmdinput = open(vmdinput_file, "w")

    Base.write(vmdinput, "mol new \"$xyzfile\" \n")
    units = split(selection, " ")
    for u in units
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $u\" ] \n")
        Base.write(vmdinput, "\$sel writepdb \"$(rawname * "_" * u * ".pdb")\" \n")
    end
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)
    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    return rawname, vmdoutput

end

function _convert_XYZ2PDB(xyzfile::String, vmdselection::String; vmd="vmd", pdbfile=nothing)
    
    if isnothing(pdbfile); pdbfile = tempname() * "pdb"; end
    vmdinput_file = tempname() * ".tcl"
    vmdinput = open(vmdinput_file, "w")

    Base.write(vmdinput, "mol new \"$xyzfile\" \n")
    Base.write(vmdinput, "set sel [ atomselect top \"$vmdselection\" ] \n")
    Base.write(vmdinput, "\$sel writepdb \"$pdbfile\" \n")
    Base.write(vmdinput, "exit \n")
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")
    
    return pdbfile, vmdoutput
end


function Z_propagation(atoms, x, y, z, zsize, phase)

    xcoords = Float64[]; ycoords = Float64[]; zcoords = Float64[]; atomnames = String[];
    parameters = get_crystallographic_info(phase)[3]
 
    # Iβ
    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
        c = parameters[1][3];
        for k in collect(1:1:zsize)
            append!(atomnames, atoms); append!(xcoords, x); append!(ycoords, y); append!(zcoords, z .+ c*(k-1));
        end        
    end

    return atomnames, xcoords, ycoords, zcoords

end


function XY_coord_trimming(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}, cell_dim::Vector{Int64}, phase::String)

    atomnames, xcoords, ycoords, zcoords = String[], Float64[], Float64[], Float64[]
    xsize, ysize = cell_dim[1], cell_dim[2]

    units = 1

    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
        for j in collect(1:1:ysize), i in collect(1:1:xsize)
            atomstemp, xtemp, ytemp, ztemp = atoms[units], x[units], y[units], z[units]
            _atomselect_indexes = (eachindex(atomstemp) .== 0)
            if ((i == 1) && (j != ysize)) || ((i != 1) && (i != xsize) && (j == 1))
                _atomselect_indexes = (eachindex(atomstemp) .>= 64) .& (eachindex(atomstemp) .<= 84)
            end
            if ((j != 1) && (i == xsize)) || ((i != 1) && (i != xsize) && (j == ysize))
                _atomselect_indexes = (eachindex(atomstemp) .>= 22) .& (eachindex(atomstemp) .<= 42)
            end
            if ((i == 1) && (j == ysize)) || ((j == 1) && (i == xsize))
                _atomselect_indexes = ((eachindex(atomstemp) .>= 22) .& (eachindex(atomstemp) .<= 42)) .| ((eachindex(atomstemp) .>= 64) .& (eachindex(atomstemp) .<= 84))
            end
            append!(atomnames, atomstemp[.!(_atomselect_indexes)])
            append!(xcoords, xtemp[.!(_atomselect_indexes)]); append!(ycoords, ytemp[.!(_atomselect_indexes)]); append!(zcoords, ztemp[.!(_atomselect_indexes)]);

            units += 1
        end
    else error("The phase $phase is not implemented yet."); end
    
    return atomnames, xcoords, ycoords, zcoords

end
