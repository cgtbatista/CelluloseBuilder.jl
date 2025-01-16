function pbcXYZ(
            xyzname::String,
            xyzsizes::Vector{Int64};
            phase="Iβ", structure=nothing, fibril=nothing, new_xyzname=joinpath(tempdir(), "cellulose.xyz"),
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
        selection = if haskey(microfibril, key)
                microfibril[key]
            else
                error("The fibril $fibril is not implemented yet on the microfibril dictionary.")
        end
    end
   
    if !isfibril
        selection = monolayer(xyzsizes, phase=phase, structure=structure)
    end

    nfragments = length(split(selection, " "))

    tcl = tempname() * ".tcl"

    vmdinput = open(tcl, "w")
    Base.write(vmdinput, """
    mol new $xyzname
    set sel [ atomselect top \"fragment $selection\" ]
    \$sel writexyz $new_xyzname
    exit
    """)
    Base.close(vmdinput)

    vmdoutput = split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if vmdDebug
        return vmdoutput
    else
        return new_xyzname, selection, nfragments
    end
end

"""
    pbcXYZ(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")

Applies the periodic boundary conditions on the cellulose crystal over *xy-plane* using the number of fragments (e.g. chains) as parameter.
It returns cellulose **XYZ file** needed to prepare the PDB and PSF files.
"""
function pbcXYZ(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")
    
    xyz = joinpath(tempdir(), "cellulose.xyz") ## the new and last XYZ file!

    a=nfrag; boundary=1; n_forbbiden=0; upper=nfrag-1;
    forbbiden=Int64[]; remainder=Int64[];

    if phase == "Ib" || phase == "Iβ" || phase == "II"

        if pbc == :all || pbc == :ALL || pbc == :All
            for b in collect(boundary:1:ysize)
                n_forbbiden = (2*xsize-1)*b - 1
                push!(forbbiden, convert(Int64, n_forbbiden))
            end
            for b in collect(boundary:1:xsize)
                a = a - 1
                push!(forbbiden, convert(Int64, a))
            end
        elseif pbc == :a || pbc == :A
            for b in collect(boundary:1:ysize)
                n_forbbiden = (2*xsize-1)*b - 1
                push!(forbbiden, convert(Int64, n_forbbiden))
            end
        elseif pbc == :b || pbc == :B
            for b in collect(boundary:1:xsize)
                a = a - 1
                push!(forbbiden, convert(Int64, a))
            end
        end
       
        for num in 0:upper
            dummy_logical = 1
            for ith_forbs in forbbiden
                if num == ith_forbs
                    dummy_logical = 0
                    break
                end
            end
            if dummy_logical == 1
                remainder = push!(remainder, convert(Int64, num))
            end
        end

        sel_fragments = join(remainder, " "); n_fragments = length(remainder);
        
        vmdinput_file = tempname() * ".tcl"
        vmdinput = open(vmdinput_file, "w")
        Base.write(vmdinput, "mol new \"$xyzfile\" \n")
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $sel_fragments\" ] \n")
        Base.write(vmdinput, "\$sel writexyz \"$xyz\" \n")
        Base.write(vmdinput, "exit \n")
        Base.close(vmdinput)
        
        vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    elseif phase == "Ia" || phase == "Iα"

        remainder = collect(0:1:(nfrag-1))
        sel_fragments = join(remainder, " "); n_fragments = length(remainder);
        
        vmdinput_file = tempname() * ".tcl"
        vmdinput = open(vmdinput_file, "w")
        Base.write(vmdinput, "mol new \"$xyzfile\" \n")
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $sel_fragments\" ] \n")
        Base.write(vmdinput, "\$sel writexyz \"$xyz\" \n")
        Base.write(vmdinput, "exit \n")
        Base.close(vmdinput)

        vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        
        remainder = collect(0:1:(nfrag-1))
        sel_fragments = join(remainder, " "); n_fragments = length(remainder);
        
        vmdinput_file = tempname() * ".tcl"
        vmdinput = open(vmdinput_file, "w")
        Base.write(vmdinput, "mol new \"$xyzfile\" \n")
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $sel_fragments\" ] \n")
        Base.write(vmdinput, "\$sel writexyz \"$xyz\" \n")
        Base.write(vmdinput, "exit \n")
        Base.close(vmdinput)
        
        vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    else
        error("The phase $phase is not implemented yet.")
    end


    return xyz, sel_fragments, n_fragments, vmdoutput

end