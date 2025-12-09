"""
    PBC(xsize::Int64, ysize::Int64, zsize::Int64; phase="Iβ", pbc=nothing)

Returns `(ncells, lattice)` where `ncells` is the number of unit cells along x, y and z axes and `lattice` is the unit cell vector of this cell.
This function is used to define the periodic boundary conditions required to build the cellulose crystal.
"""
function PBC(
    xsize::T, ysize::T, zsize::T; phase="Iβ", pbc=nothing
) where T
    phase = lowercase(phase)
    pbc = isnothing(pbc) ? nothing : lowercase(pbc)
    valid_phases = in(
        phase,
        Set(["iβ", "ib", "iα", "ia", "ii", "iii", "iii_i", "iiii"])
    )
    valid_pbcs = in(
        pbc,
        Set(["a", "b", "all"])
    ) || isnothing(pbc)
    if !(valid_phases || valid_pbcs)
        throw(ArgumentError("The phase $phase or the pbc $pbc is not valid."))
    end
    isalpha = in(phase, ["iα", "ia"])
    a, b, c = isalpha ? (zsize, ysize, xsize) : (xsize, ysize, zsize)
    if isnothing(pbc)
        println("""
        Periodic boundary conditions will not be special applied.
        Default translational symmetry will be applied with the surfaces (1 0 0) and (0 1 0) exposed!
        """)
    end
    if !isnothing(pbc)
        dictionary = pbcdata(xsize, ysize, zsize)
        key = (phase, lowercase(pbc))
        if !haskey(dictionary, key)
            throw(ArgumentError("The `$pbc` periodic boundary condition is not possible for the $phase phase."))
        end
        xsize, ysize, zsize, pbc_message = dictionary[key]()
        print(pbc_message)
    end
    return [xsize, ysize, zsize], [a, b, c]
end

"""
    PBC(xyzsizes::Vector{Int64}; phase="Iβ", pbc=nothing)

See `PBC(xsize::Int64, ysize::Int64, zsize::Int64; phase="Iβ", pbc=nothing)`.
"""
function PBC(xyzsizes::Vector{T}; phase="Iβ", pbc=nothing) where T
    return PBC(xyzsizes[1], xyzsizes[2], xyzsizes[3], phase=phase, pbc=pbc)
end

"""
    PBC(monomers::Int64; phase="Iβ")

Define the crystallographic space to build the *cellulose fibril crystal* with a degree of β-glucose units (`monomers`).
"""
function PBC(monomers::Int64; phase="Iβ")
    n = Int64(monomers/2)
    isalpha = in(lowercase(phase), ["iα", "ia"])
    a, b, c = isalpha ? (n, 0, 0) : (0, 0, n)
    dictionary, key = pbcdata(monomers), lowercase(phase)

    xsize, ysize, zsize = if haskey(dictionary, key)
             dictionary[key]
        else
            error("The $phase phase is not possible to build.")
    end

    println("""
    Periodic boundary conditions will be applied to build the $phase elementary fibril.
    There is $monomers of β-glucose monomers on each cellulose chain — $n cellobiose units per chain.
    """)

    return [xsize, ysize, zsize], [a, b, c]
end

"""
    PBC(monomers::Int64, chains::Int64; phase="Iβ", layer="center")

Define the crystallographic space to craft a *cellulose layer* with a degree of β-glucose units (`monomers`) and a number of cellulose chains (`chains`).
"""
function PBC(monomers::Int64, chains::Int64; phase="Iβ", layer="center")

    n = Int64(monomers/2)

    a, b, c = if lowercase(phase) in ["ib", "iβ"]
        1, chains, n
    elseif lowercase(phase) in ["ia", "iα"]
        n, 0, 0
    elseif lowercase(phase) == "ii"
        chains, 1, n
    else
        error("The $phase phase is not available or possible to use in this building.")
    end

    dictionary, key = pbcdata(monomers, chains=chains), (lowercase(phase), lowercase(layer))
    
    xsize, ysize, zsize = if haskey(dictionary, key)
            dictionary[key]
        else
            error("The layer `$layer` is not possible to build for $phase phase.")
    end
    
    println("""
    Periodic boundary conditions will be applied to build the `$(lowercase(layer))` layer from $phase polymorph.
    Building monolayer composed of $chains cellulose chains with $n cellobiose units ($monomers monomers).
    """)

    return [xsize, ysize, zsize], [a, b, c]
end

"""
    pbcdata(xsize::Int64, ysize::Int64, zsize::Int64)

Dictionary useful for `(a, b, c)` definition of the crystal.
"""
function pbcdata(xsize::Int64, ysize::Int64, zsize::Int64)
    
    PBC = Dict{Tuple{String, String}, Function}(
        ("ib", "a") => () -> (xsize + 1, ysize, zsize, """
        Periodic boundary conditions will be applied on crystalographic direction `a`.
        Surfaces (1 0 0), (2 0 0), and (0 1 0) will be exposed!
        
        """),
        ("ib", "b") => () -> (xsize, ysize + 1, zsize, """
        Periodic boundary conditions will be applied on crystalographic direction `b`.
        Surfaces (1 0 0), (0 1 0), and (0 2 0) will be exposed!
        
        """),
        ("ib", "all") => () -> (xsize + 1, ysize + 1, zsize, """
        Periodic boundary conditions will be applied on both crystalographic directions `a` and `b`.
        Surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!
        
        """),
        ("ii", "a") => () -> (xsize + 1, ysize, zsize, """
        Periodic boundary conditions will be applied on crystalographic direction `a`.
        Surfaces (1 0 0), (2 0 0), and (0 1 0) will be exposed!
        
        """),
        ("ii", "b") => () -> (xsize, ysize + 1, zsize, """
        Periodic boundary conditions will be applied on crystalographic direction `b`.
        Surfaces (1 0 0), (0 1 0), and (0 2 0) will be exposed!
        
        """),
        ("ii", "all") => () -> (xsize + 1, ysize + 1, zsize, """
        Periodic boundary conditions will be applied on both crystalographic directions `a` and `b`.
        Surfaces (1 0 0), (2 0 0), (0 1 0), and (0 2 0) will be exposed!
        
        """)
    )

    aliases = [
        (("iβ", "a"), ("ib", "a")),
        (("iβ", "b"), ("ib", "b")),
        (("iβ", "all"), ("ib", "all"))
    ]

    for (alias, original) in aliases
        PBC[alias] = PBC[original]
    end

    return PBC
end

"""
    pbcdata(monomers::Int64; chains=nothing)

Dictionary useful for cellulose definition of the layers and fibril.
"""
function pbcdata(monomers::Int64; chains=nothing)
    
    n = Int64(monomers/2)

    if isnothing(chains)
        PBC = Dict{String, Tuple{Int64, Int64, Int64}}(
            "ib" => (5, 7, n),
            "II" => (7, 5, n),
            "ia" => (7, 6, n)
        )

        aliases = [
            ("iβ", "ib"),
            ("iα", "ia")
        ]
    end

    if typeof(chains) == Int64
        chains2 = chains + 1        
        PBC = Dict{Tuple{String, String}, Tuple{Int64, Int64, Int64}}(
            ("ib", "center") => (2, chains2, n),
            ("ib", "origin") => (2, chains, n),
            ("ia", "monolayer") => (chains, chains, n),
            ("ii", "center") => (chains2, 2, n),
            ("ii", "origin") => (chains, 2, n)
        )

        aliases = [
            (("ib", "c"), ("ib", "center")),
            (("ib", "o"), ("ib", "origin")),
            (("ia", "m"), ("ia", "monolayer")),
            (("ii", "c"), ("ii", "center")),
            (("ii", "o"), ("ii", "origin")),
            (("iβ", "center"), ("ib", "center")),
            (("iβ", "origin"), ("ib", "origin")),
            (("iα", "monolayer"), ("ia", "monolayer")),
            (("iβ", "c"), ("ib", "center")),
            (("iβ", "o"), ("ib", "origin")),
            (("iα", "m"), ("ia", "monolayer"))
        ]
    end

    for (alias, original) in aliases
        PBC[alias] = PBC[original]
    end

    return PBC
end