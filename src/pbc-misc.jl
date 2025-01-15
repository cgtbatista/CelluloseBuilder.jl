"""
    pbcdata(xsize::Int64, ysize::Int64, zsize::Int64)

Dictionary useful for ```(a,b,c)``` definition of the crystal.
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

const microfibril = Dict{Tuple{String, String}, String}([
        [("iβ", "34566543"), ("ib", "34566543")] .=> "18 27 36 10 19 28 37 11 20 29 38 47 3 12 21 30 39 48 4 13 22 31 40 49 5 14 23 32 41 15 24 33 42 16 25 34";
        [("iβ", "33333333"), ("ib", "33333333")] .=> "18 27 36 10 19 28 20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34";
        [("iβ", "23454321"), ("ib", "23454321")] .=> "27 36 19 28 37 20 29 38 47 12 21 30 39 48 22 31 40 49 23 32 41 33 42 34";
        [("iβ", "12333321"), ("ib", "12333321")] .=> "18 10 19 11 20 29 12 21 30 22 31 40 23 32 41 33 42 34";
        [("iβ", "234432"), ("ib", "234432")]     .=> "19 28 20 29 38 12 21 30 39 13 22 31 40 14 23 32 24 33";
        [("iβ", "333333"), ("ib", "333333")]     .=> "20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34";
        [("iα", "34566543"), ("ia", "34566543")] .=> "2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 38 39";
        [("iα", "33333333"), ("ia", "33333333")] .=> "2 7 8 9 12 13 14 15 16 18 19 20 21 22 23 25 26 27 28 29 32 33 34 39";
        [("iα", "23454321"), ("ia", "23454321")] .=> "12 13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33 34 36 37 38 39";
        [("iα", "12333321"), ("ia", "12333321")] .=> "7 8 9 13 14 15 19 20 21 25 26 27 31 32 33 37 38 39";
        [("iα", "234432"), ("ia", "234432")]     .=> "13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33";
        [("iα", "333333"), ("ia", "333333")]     .=> "14 15 19 20 21 22 24 25 26 27 28 29 31 32 33 34 38 39";
        [("ii", "34566543")]                     .=> "3 4 5 6 7 8 9 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 41 42 43 44 45 46 47";
])

function monolayer(xyzsizes::Vector{Int64}; phase="Iα", structure="monolayer")

    selection = ""

    if in(lowercase(structure), Set(["monolayer", "m"])) && in(lowercase(phase), Set(["iα", "ia"]))
        nsize = xyzsizes[2]
        for chain in 1:nsize
            selection = selection * " " * string(chain*(nsize-1))
        end

        return strip(selection)
    end

    isorigin, iscenter = in(lowercase(structure), Set(["origin", "o"])), in(lowercase(structure), Set(["center", "c"]))
    
    if in(lowercase(phase), Set(["ib", "iβ"]))
        nsize = xyzsizes[2]
        if isorigin
            for chain in 1:nsize
                selection = selection * " " * string(3*(chain-1))
            end
        end
        if iscenter
            for chain in 1:nsize
                if chain > nsize-1; continue; end
                selection = selection * " " * string(3*(chain-1)+1)
            end
        end

        return strip(selection)
    end

    if lowercase(phase) == "ii"
        nsize = xyzsizes[1]
        if isorigin
            for chain in 1:nsize
                selection = selection * " " * string(2*(chain-1))
            end
        end
        if iscenter
            for chain in 1:nsize
                if chain > nsize-1; continue; end
                selection = selection * " " * string(2*(chain-1)+1)
            end
        end

        return strip(selection)
    end

    if selection == ""
        error("The phase $phase does not supports the structure $structure.")
    end
end

function transformingPBC(
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
    transformingPBC(nfrag::Int64, xsize::Int64, ysize::Int64; phase="Iβ", pbc=nothing, xyzfile="file.xyz", vmd="vmd")

Applies the periodic boundary conditions on the cellulose crystal over *xy-plane* using the number of fragments (e.g. chains) as parameter.
It returns cellulose **XYZ file** needed to prepare the PDB and PSF files.
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