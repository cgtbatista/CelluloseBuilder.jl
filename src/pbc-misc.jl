function _PBC_conditional_center(xyz::String, nsize::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd", phase="Iβ")

    for chain in collect(1:1:nsize)
        if chain > nsize-1; continue; end
        if phase == "Ib" || phase == "Iβ"
            n = 3*(chain-1)+1
        elseif phase == "II"
            n = 2*(chain-1)+1
        end
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

function _PBC_conditional_origin(xyz::String, nsize::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd", phase="Iβ")

    for chain in collect(1:1:nsize)
        if phase == "Ib" || phase == "Iβ"
            n = 3*(chain-1)
        elseif phase == "II"
            n = 2*(chain-1)
        end
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

function _PBC_conditional_monolayer(xyz::String, nsize::Int64, remainder::Vector{Int64}; xyzfile="filename.xyz", vmd="vmd")

    for chain in collect(1:1:nsize)
        n = chain*(nsize-1)
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

function _PBC_conditional_fibril(xyz::String; xyzfile="filename.xyz", vmd="vmd", phase="Iβ", fibril=nothing)

    if phase == "Iβ" || phase == "Ib"
        if isnothing(fibril)
            remainder = split("18 27 36 10 19 28 37 11 20 29 38 47 3 12 21 30 39 48 4 13 22 31 40 49 5 14 23 32 41 15 24 33 42 16 25 34", " ")
        elseif fibril == "234432"
            remainder = split("19 28 20 29 38 12 21 30 39 13 22 31 40 14 23 32 24 33", " ")
        elseif fibril == "333333"
            remainder = split("20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34", " ")
        elseif fibril == "12333321"
            remainder = split("18 10 19 11 20 29 12 21 30 22 31 40 23 32 41 33 42 34", " ")
        elseif fibril == "33333333"
            remainder = split("18 27 36 10 19 28 20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34", " ")
        elseif fibril == "23454321"
            remainder = split("27 36 19 28 37 20 29 38 47 12 21 30 39 48 22 31 40 49 23 32 41 33 42 34", " ")
        end
    elseif phase == "II"
        if isnothing(fibril)
            remainder = split("3 4 5 6 7 8 9 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 41 42 43 44 45 46 47", " ")
        end
    elseif phase == "Iα" || phase == "Ia"
        if isnothing(fibril)
            remainder = split("2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 38 39", " ")
        elseif fibril == "234432"
            remainder = split("13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33", " ")
        elseif fibril == "333333"
            remainder = split("14 15 19 20 21 22 24 25 26 27 28 29 31 32 33 34 38 39", " ")
        elseif fibril == "12333321"
            remainder = split("7 8 9 13 14 15 19 20 21 25 26 27 31 32 33 37 38 39", " ")
        elseif fibril == "33333333"
            remainder = split("2 7 8 9 12 13 14 15 16 18 19 20 21 22 23 25 26 27 28 29 32 33 34 39", " ")
        elseif fibril == "23454321"
            remainder = split("12 13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33 34 36 37 38 39", " ")
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

function transformingPBC(style::String, xyzsize::Vector{Int64}; phase="Iβ", fibril=nothing, xyzfile="filename.xyz", vmd="vmd", vmdDebug=false)
    
    xyz = joinpath(tempdir(), "cellulose.xyz") ## the new and last XYZ file!

    boundary=1; remainder=Int64[];
    
    if phase == "Ib" || phase == "Iβ"

        if style == "center"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_center(xyz, xyzsize[2], remainder, xyzfile=xyzfile, vmd=vmd, phase=phase)
        elseif style == "origin"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_origin(xyz, xyzsize[2], remainder, xyzfile=xyzfile, vmd=vmd, phase=phase)
        elseif style == "fibril"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_fibril(xyz, xyzfile=xyzfile, vmd=vmd, fibril=fibril, phase=phase)
        else
            error("The phase $phase does not supports $style, only the styles center, origin and fibril.")
        end

    elseif phase == "II"

        if style == "center"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_center(xyz, xyzsize[1], remainder, xyzfile=xyzfile, vmd=vmd, phase=phase)
        elseif style == "origin"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_origin(xyz, xyzsize[1], remainder, xyzfile=xyzfile, vmd=vmd, phase=phase)
        elseif style == "fibril"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_fibril(xyz, xyzfile=xyzfile, vmd=vmd, fibril=fibril, phase=phase)
        else
            error("The phase $phase does not supports $style, only the styles center, origin and fibril.")
        end

    elseif phase == "Ia" || phase == "Iα"

        if style == "monolayer"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_monolayer(xyz, xyzsize[2], remainder, xyzfile=xyzfile, vmd=vmd)
        elseif style == "fibril"
            xyz, sel_fragments, n_fragments, vmdoutput = _PBC_conditional_fibril(xyz, xyzfile=xyzfile, vmd=vmd, fibril=fibril, phase=phase)
        else
            error("The phase $phase does not supports $style, only the styles monolayer and fibril.")
        end

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        
        error("The phase $phase does not supports $style or any other style.")

    else
        error("The phase $phase is not implemented yet.")
    end

    if vmdDebug
        return vmdoutput
    else
        return xyz, sel_fragments, n_fragments
    end
end
