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

    if isnothing(fibril)
        if phase == "Iβ" || phase == "Ib"
            remainder = "18 27 36 10 19 28 37 11 20 29 38 47 3 12 21 30 39 48 4 13 22 31 40 49 5 14 23 32 41 15 24 33 42 16 25 34"
        elseif phase == "II"
            remainder = "3 4 5 6 7 8 9 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 41 42 43 44 45 46 47"
        elseif phase == "Iα" || phase == "Ia"
            remainder = "2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 38 39"
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

function transformingPBC(style::String, xyzsize::Vector{Int64}; phase="Iβ", fibril=nothing, xyzfile="filename.xyz", vmd="vmd")
    
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


    return xyz, sel_fragments, n_fragments, vmdoutput

end


"""

    gettingPBC(xsize::Int64, ysize::Int64, zsize::Int64, phase::String; pbc=nothing)
    gettingPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    gettingPBC(style::String, j::Int64, k::Int64, phase::String)
    gettingPBC(units::Int64, phase::String)

This function aims to return the crystallographic information for the setted cellulose phase. The information is related to the CHARMM atomnames, the unit cell
parameters and the fractional coordinates of the asymetric unit.

## Arguments

- `xyzsizes::Vector{Int64}`: The number of unit cells along x, y and z axes (`a`, `b` and `c`).
- `xsize::Int64`: The number of unit cells along x axis (`a`).
- `ysize::Int64`: The number of unit cells units along y axis (`b`).
- `zsize::Int64`: The number of unit cells units along z axis (`c`).
- `style::String`: Type of special cellulose monolayer (monolayer, center or origin).
- `units::Int64`: The number of fibril cellobiose units.
- `j::Int64` and `k::Int64`: The number of unit cell units.
- `phase::String`: The cellulose phase. It could be `Iβ`, `Iα`, `II` or `III`.
- `pbc=nothing`: The periodic boundary conditions to be applied. It's could be around `:A`, `:B`, or both directions `:ALL`. The default is `nothing`.

### Examples

```jldoctest

julia > gettingPBC(5, 7, 8, "Iβ")
julia > gettingPBC([ 5, 7, 8 ], "Iβ")

```

"""


function gettingPBC(xsize::Int64, ysize::Int64, zsize::Int64, phase::String; pbc=nothing)
    if phase == "Ib" || phase == "Iβ"
        if pbc == :all || pbc == :All || pbc == :ALL
            xsize += 1; ysize += 1;
            println("         periodic boundary conditions will be applied in a and b crystalographic directions.")
            println("         surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!")
        elseif pbc == :a || pbc == :A
            xsize += 1;
            println("         periodic boundary conditions will be applied in a crystalographic direction.")
            println("         surfaces (1 0 0), (2 0 0), and (0 1 0) will be exposed!")
        elseif pbc == :b || pbc == :B
            ysize += 1;
            println("         periodic boundary conditions will be applied in b crystalographic direction.")
            println("         surfaces (1 0 0), (0 1 0), and (0 2 0) will be exposed!")
        elseif isnothing(pbc)
            println("         periodic boundary conditions will not be special applied.")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        end
    elseif phase == "Ia" || phase == "Iα"
        if !isnothing(pbc)
            println("         periodic boundary conditions $pbc will not be applied in the Iα phase, because it is not valid!")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        else
            println("         periodic boundary conditions will not be special applied.")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        end
    elseif phase == "II"
        if pbc == :all || pbc == :All || pbc == :ALL
            xsize += 1;
            println("         periodic boundary conditions will be applied in a and b crystalographic directions.")
            println("         surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!")
        elseif pbc == :a || pbc == :A
            xsize += 1;
            println("         periodic boundary conditions will be applied in a crystalographic direction.")
            println("         surfaces (1 0 0), (2 0 0), and (0 1 0) will be exposed!")
        elseif pbc == :b || pbc == :B
            ysize += 1;
            println("         periodic boundary conditions will be applied in b crystalographic direction.")
            println("         surfaces (1 0 0), (0 1 0), and (0 2 0) will be exposed!")
        elseif isnothing(pbc)
            println("         periodic boundary conditions will not be special applied.")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        end
    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        if !isnothing(pbc)
            println("         periodic boundary conditions $pbc will not be applied in the Iα phase, because it is not valid!")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        else
            println("         periodic boundary conditions will not be special applied.")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        end
    else
        error("The phase $phase is not implemented yet.")
    end
    return [xsize, ysize, zsize]
end

function gettingPBC(units::Int64, phase::String)
    
    if phase == "Ib" || phase == "Iβ" || phase == "II"

        xsize = 5; ysize = 7; zsize = units;
        xlattice = 0; ylattice = 0; zlattice = zsize;
        
    elseif phase == "Ia" || phase == "Iα"

        xsize = 7; ysize = 6; zsize = units;
        xlattice = zsize; ylattice = 0; zlattice = 0;

    else
        error("The phase $phase cannot be used as a base for elementary fibril.")
    end

    println("         periodic boundary conditions will be applied to build the elementary fibril.")
    println("         number of cellobiose residues per cellulose chain is $zsize.")

    return [xsize, ysize, zsize], [xlattice, ylattice, zlattice]
end


function gettingPBC(style::String, j::Int64, k::Int64, phase::String)
    
    if phase == "Ib" || phase == "Iβ"

        if style == "center"
            xsize = 2; ysize = j+1; zsize = k;
            xlattice = 1; ylattice = j; zlattice = k;
            println("         building monolayer composed of $j cellulose units.")
            println("         `center` labelled cellulose with $zsize cellobiose units!")
        elseif style == "origin"
            xsize = 2; ysize = j; zsize = k;
            xlattice = 1; ylattice = j; zlattice = k;
            println("         building monolayer composed of $j cellulose units.")
            println("         `origin` labelled cellulose with $zsize cellobiose units!")
        else
            error("The option $style is not available for the $phase phase or is not available at all...")
        end

    elseif phase == "Ia" || phase == "Iα"
        
        if style == "monolayer"
            xsize = j; ysize = j; zsize = k;
            xlattice = zsize; ylattice = 0; zlattice = 0;
            println("         building monolayer composed of $j cellulose units.")
            println("         `monolayer` labelled cellulose with $zsize cellobiose units!")
        else 
            error("The option $style is not available for the $phase phase or is not available at all...")
        end

    elseif phase == "II"

        if style == "center"
            xsize = j+1; ysize = 2; zsize = k;
            xlattice = j; ylattice = 1; zlattice = k;
            println("         building monolayer composed of $j cellulose units.")
            println("         `center` labelled cellulose with $zsize cellobiose units!")
        elseif style == "origin"
            xsize = j; ysize = 2; zsize = k;
            xlattice = j; ylattice = 1; zlattice = k;
            println("         building monolayer composed of $j cellulose units.")
            println("         `origin` labelled cellulose with $zsize cellobiose units!")
        else 
            error("The option $style is not available for the $phase phase or is not available at all...")
        end

    else
        error("The style $style is not implemented on $phase cellulose polymorph.")
    end
    return [xsize, ysize, zsize], [xlattice, ylattice, zlattice]
end

function gettingPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    xsize = xyzsizes[1]; ysize = xyzsizes[2]; zsize = xyzsizes[3];
    return gettingPBC(xsize, ysize, zsize, phase; pbc=pbc)
end