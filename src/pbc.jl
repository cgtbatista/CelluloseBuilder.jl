"""

    gettingPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    gettingPBC(xsize::Int64, ysize::Int64, zsize::Int64, phase::String; pbc=nothing)

This function aims to return the crystallographic information for the setted cellulose phase. The information is related to the CHARMM atomnames, the unit cell
parameters and the fractional coordinates of the asymetric unit.

## Arguments

- `xyzsizes::Vector{Int64}`: The number of unit cells along x, y and z axes (`a`, `b` and `c`).
- `xsize::Int64`: The number of unit cells along x axis (`a`).
- `ysize::Int64`: The number of unit cells units along y axis (`b`).
- `zsize::Int64`: The number of unit cells units along z axis (`c`).
- `phase::String`: The cellulose phase. It could be `Iβ`, `Iα`, `II` or `III`.
- `pbc=nothing`: The periodic boundary conditions to be applied. It's could be around `:A`, `:B`, or both directions `:ALL`. The default is `nothing`.

### Examples

```julia-repl

julia > 

```

"""

function gettingPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    xsize = xyzsizes[1]; ysize = xyzsizes[2]; zsize = xyzsizes[3];
    return gettingPBC(xsize, ysize, zsize, phase; pbc=pbc)
end

function gettingPBC(xsize::Int64, ysize::Int64, zsize::Int64,  phase::String; pbc=nothing)
    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
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
    elseif phase == "I-ALPHA" || phase == "Ia" || phase == "Iα"
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
    elseif phase == "III"
        if !isnothing(pbc)
            println("         periodic boundary conditions $pbc will not be applied in the Iα phase, because it is not valid!")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        else
            println("         periodic boundary conditions will not be special applied.")
            println("         default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!")
        end
    end
    return [xsize, ysize, zsize]
end

function PBCtools(phase, pbc, nfrag, xsize, ysize, xyzfile, vmd)
    
    forbbiden=Int64[]; remainder=Int64[];

    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"

        a=nfrag; boundary=1; n_forbbiden=0; upper=nfrag-1;
        
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

        new_xyzfile = "/tmp/cellulose" * ".xyz"
        vmdinput_file = tempname() * ".tcl"
        
        vmdinput = open(vmdinput_file, "w")
        Base.write(vmdinput, "mol new \"$xyzfile\" \n")
        Base.write(vmdinput, "set sel [ atomselect top \"fragment $sel_fragments\" ] \n")
        Base.write(vmdinput, "\$sel writexyz \"$new_xyzfile\" \n")
        Base.write(vmdinput, "exit \n")
        Base.close(vmdinput)
        vmdoutput = split(Base.read(`$vmd -dispdev text -e $vmdinput_file`, String), "\n")

    end

    return new_xyzfile, sel_fragments, n_fragments

end