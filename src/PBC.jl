function _PBC_conditional_center(xyz::String, boundary::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd", nsize=1)

    for b in collect(boundary:1:nsize)
        if b > nsize-1; continue; end
        n = 2*(b-1)+1
        push!(remainder, convert(Int64, n))
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

    return xyz, sel_fragments, n_fragments, vmdoutput
end

function _PBC_conditional_origin(xyz::String, boundary::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd", nsize=1)

    for b in collect(boundary:1:nsize)
        n = 2*(b-1)
        push!(remainder, convert(Int64, n))
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

    return xyz, sel_fragments, n_fragments, vmdoutput
end

function _PBC_conditional_monolayer(xyz::String, boundary::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd", nsize=1)

    for b in collect(boundary:1:nsize)
        n = b*(nsize-1)
        push!(remainder, convert(Int64, n))
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

    return xyz, sel_fragments, n_fragments, vmdoutput
end

function _PBC_conditional_fibril(xyz::String; xyzfile="filename.xyz", vmd="vmd", fibril=nothing)

    if isnothing(fibril)
        remainder = "3 4 5 6 7 8 9 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 41 42 43 44 45 46 47"
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

    return xyz, sel_fragments, n_fragments, vmdoutput
end


"""

    transformingPBC(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")
    transformingPBC(style::String, nsize::Int64; phase="Iβ", fibril=nothing, xyzfile="file.xyz", vmd="vmd")

This function applies the periodic boundary conditions on the cellulose fibrils over the axis `a`, `b`, and `both`.
After this transformation, you export the cellulose fibril XYZ file with the cartesian coordinates.

## Arguments

- `nfrag::Int64`: is the number of fragments inside the XYZ file.
- `xsize::Int64`: the number of unit cells along x axis (`a`).
- `ysize::Int64`: the number of unit cells units along y axis (`b`).

### Examples

```jldoctest

julia > transformingPBC(59, 5, 7)

```

"""

function transformingPBC(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")
    
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

        mv(xyzfile, xyz, force=true)

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        
        mv(xyzfile, xyz, force=true)

    else
        error("The phase $phase is not implemented yet.")
    end


    return xyz, sel_fragments, n_fragments, vmdoutput

end

function transformingPBC(style::String, nsize::Int64; phase="Iβ", fibril=nothing, xyzfile="file.xyz", vmd="vmd")
    
    xyz = joinpath(tempdir(), "cellulose.xyz") ## the new and last XYZ file!

    boundary=1; remainder=Int64[];
    
    if phase == "Ib" || phase == "Iβ" || phase == "II"

        if style == "center"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_center(xyz, boundary, remainder, xyzfile=xyzfile, vmd=vmd, nsize=nsize)
        elseif sytle == "origin"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_origin(xyz, boundary, remainder, xyzfile=xyzfile, vmd=vmd, nsize=nsize)
        elseif style == "fibril"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_fibril(xyz, xyzfile=xyzfile, vmd=vmd, fibril=fibril)
        else
            error("The phase $phase does not supports $style, only the styles center, origin and fibril.")
        end

    elseif phase == "Ia" || phase == "Iα"

        if style == "monolayer"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_monolayer(xyz, boundary, remainder, xyzfile=xyzfile, vmd=vmd, nsize=nsize)
        elseif style == "fibril"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_fibril(xyz, xyzfile=xyzfile, vmd=vmd, fibril=fibril)
        else
            error("The phase $phase does not supports $style, only the styles monolayer and fibril.")
        end

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        
        error("The phase $phase does not supports $style or any other style.")

    else
        error("The phase $phase is not implemented yet.")
    end


    return xyz, sel_fragments, n_fragments, vmdoutput

end