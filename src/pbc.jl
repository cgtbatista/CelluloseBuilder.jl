"""

    getPBC(xsize::Int64, ysize::Int64, zsize::Int64, phase::String; pbc=nothing)
    getPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    getPBC(style::String, j::Int64, k::Int64, phase::String)
    getPBC(units::Int64, phase::String)

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

julia > getPBC(5, 7, 8, "Iβ")
julia > getPBC([ 5, 7, 8 ], "Iβ")

```

"""
function getPBC(xsize::Int64, ysize::Int64, zsize::Int64, phase::String; pbc=nothing)

    valid_phases = Set(["iβ", "ib", "iα", "ia", "ii", "iii", "iii_i", "iiii"])

    if !in(lowercase(phase), valid_phases)
        error("The phase $phase is not implemented yet.")
    end

    a, b, c = if lowercase(phase) in ["iα", "ia"]
        copy(zsize), copy(ysize), copy(xsize)
    else
        copy(xsize), copy(ysize), copy(zsize)
    end

    valid_pbc = Set(["a", "b", "all"])
    if !isnothing(pbc) && !in(lowercase(pbc), valid_pbc)
        error("The periodic boundary conditions $pbc is not valid.")
    end

    if isnothing(pbc)
        println("""
        Periodic boundary conditions will not be special applied.
        Default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!
        """)
    elseif lowercase(pbc) == "all"
        if lowercase(phase) in ["ib", "iβ"]
            xsize += 1; ysize += 1;
            println("""
            Periodic boundary conditions will be applied both on crystalographic directions `a` and `b`.
            surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!
            """)
        elseif lowercase(phase) == "ii"
            xsize += 1
            println("""
            Periodic boundary conditions will be applied in a and b crystalographic directions.
            Surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!
            """)
        else
            error("The `$pbc` periodic boundary condition is not possible for the $phase phase.")
        end
    elseif lowercase(pbc) == "a"
        if lowercase(phase) in ["ib", "iβ", "ii"]
            xsize += 1;
            println("""
            Periodic boundary conditions will be applied on crystalographic direction `a`.
            surfaces (1 0 0), (2 0 0), and (0 1 0) will be exposed!
            """)
        else
            error("The `$pbc` periodic boundary condition is not possible for the $phase phase.")
        end
    elseif lowercase(pbc) == "b"
        if lowercase(phase) in ["ib", "iβ", "ii"]
            ysize += 1;
            println("""
            Periodic boundary conditions will be applied on crystalographic direction `b`.
            Surfaces (1 0 0), (0 1 0), and (0 2 0) will be exposed!
            """)
        else
            error("The `$pbc` periodic boundary condition is not possible for the $phase phase.")
        end
    end

    return [xsize, ysize, zsize], [a, b, c]
end

function getPBC(xyzsizes::Vector{Int64}, phase::String; pbc=nothing)
    return getPBC(xyzsizes[1], xyzsizes[2], xyzsizes[3], phase; pbc=pbc)
end

function getPBC(units::Int64, phase::String)
    
    xsize, ysize, zsize, a, b, c = if lowercase(phase) in ["ib", "iβ"]
            5, 7, units, 0, 0, units
        elseif lowercase(phase) == "ii"
            7, 5, units, 0, 0, units
        elseif lowercase(phase) in ["ia", "iα"]
            7, 6, units, units, 0, 0
        else
            error("The phase $phase cannot be used as a base for elementary fibril.")
    end

    println("""
    Periodic boundary conditions will be applied to build the $phase elementary fibril.
    Number of cellobiose residues per cellulose chain is $zsize.
    """)

    return [xsize, ysize, zsize], [a, b, c]
end


function getPBC(layer::String, j::Int64, k::Int64, phase::String)

    xsize, ysize, zsize, a, b, c = if lowercase(phase) in ["ib", "iβ"]
            if lowercase(layer) in ["c", "center"]
                2, j+1, k, 1, j, k
            elseif lowercase(layer) in ["o", "origin"]
                2, j, k, 1, j, k
            else
                error("The option $layer is not available for the $phase phase.")
            end
        elseif lowercase(phase) in ["ia", "iα"]
            if lowercase(layer) in ["m", "monolayer"]
                j, j, k, k, 0, 0
            else
                error("The option $layer is not available for the $phase phase.")
            end
        elseif lowercase(phase) == "ii"
            if lowercase(layer) in ["c", "center"]
                j+1, 2, k, j, 1, k
            elseif lowercase(layer) in ["o", "origin"]
                j, 2, k, j, 1, k
            else
                error("The option $layer is not available for the $phase phase.")
            end
        else
            error("The phase $phase cannot be used to get this layers.")
    end
    
    println("""
    Periodic boundary conditions will be applied to build the `$(lowercase(layer))` layer.
    Building monolayer composed of $j cellulose chains with $k cellobiose units!
    """)

    return [xsize, ysize, zsize], [a, b, c]
end