const microfibril = Dict{Tuple{String, String}, String}([
        [("iβ", "34566543"), ("ib", "34566543")] .=> "18 27 36 10 19 28 37 11 20 29 38 47 3 12 21 30 39 48 4 13 22 31 40 49 5 14 23 32 41 15 24 33 42 16 25 34";
        [("iβ", "33333333"), ("ib", "33333333")] .=> "18 27 36 10 19 28 20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34";
        [("iβ", "23454321"), ("ib", "23454321")] .=> "27 36 19 28 37 20 29 38 47 12 21 30 39 48 22 31 40 49 23 32 41 33 42 34";
        [("iβ", "12333321"), ("ib", "12333321")] .=> "18 10 19 11 20 29 12 21 30 22 31 40 23 32 41 33 42 34";
        [("iβ", "234432"), ("ib", "234432")]     .=> "19 28 20 29 38 12 21 30 39 13 22 31 40 14 23 32 24 33";
        [("iβ", "23432"), ("ib", "23432")]       .=> "19 28 20 29 38 12 21 30 39 22 31 40 23 32";
        [("iβ", "3443"), ("ib", "3443")]         .=> "20 29 38 12 21 30 39 13 22 31 40 14 23 32";
        [("iβ", "333333"), ("ib", "333333")]     .=> "20 29 38 12 21 30 22 31 40 14 23 32 24 33 42 16 25 34";
        [("iα", "34566543"), ("ia", "34566543")] .=> "2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 36 37 38 39";
        [("iα", "33333333"), ("ia", "33333333")] .=> "2 7 8 9 12 13 14 15 16 18 19 20 21 22 23 25 26 27 28 29 32 33 34 39";
        [("iα", "23454321"), ("ia", "23454321")] .=> "12 13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33 34 36 37 38 39";
        [("iα", "12333321"), ("ia", "12333321")] .=> "7 8 9 13 14 15 19 20 21 25 26 27 31 32 33 37 38 39";
        [("iα", "234432"), ("ia", "234432")]     .=> "13 14 15 16 18 19 20 21 22 24 25 26 27 28 30 31 32 33";
        [("iα", "333333"), ("ia", "333333")]     .=> "14 15 19 20 21 22 24 25 26 27 28 29 31 32 33 34 38 39";
        [("ii", "34566543")]                     .=> "3 4 5 6 7 8 9 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36 41 42 43 44 45 46 47";
])

"""
    monolayer(xyzsizes::Vector{Int64}; phase="Iα", structure="monolayer")

Returns the fragments selection needed to build the monolayer structure.
"""
function monolayer(xyzsizes::Vector{Int64}; phase="Iα", structure="monolayer")

    
    if in(lowercase(structure), Set(["monolayer", "m"])) && in(lowercase(phase), Set(["iα", "ia"]))
        v_raw = 1:xyzsizes[2]
        return v_raw * (xyzsizes[2]-1)
    end

    isorigin, iscenter = in(lowercase(structure), Set(["origin", "o"])), in(lowercase(structure), Set(["center", "c"]))
    
    if in(lowercase(phase), Set(["ib", "iβ"]))
        if isorigin
            v_raw = 1:xyzsizes[2]
            return 3 * (v_raw .- 1)
        end
        if iscenter
            v_raw = 1:xyzsizes[2]-1
            return 3 * (v_raw .- 1) .+ 1
        end
    end

    if lowercase(phase) == "ii"
        if isorigin
            v_raw = 1:xyzsizes[1]
            return 2 * (v_raw .- 1)
        end
        if iscenter
            v_raw = 1:xyzsizes[1]-1
            return 2 * (v_raw .- 1) .+ 1
        end
    end
end

"""
    forbbiden(nfrags, xyzsizes::Vector{Int64}; pbc=nothing)

Returns the chains that are forbidden to be present on the crystal due to the PBCs.
"""
function forbbiden(nfrags, xyzsizes::Vector{Int64}; pbc=nothing)

    chains = Int64[]
    
    if isnothing(pbc)
        return chains
    end

    xsize, ysize = xyzsizes[1], xyzsizes[2]

    if in(lowercase(pbc), Set(["a", "all"]))
        for b in 1:ysize
            push!(chains, (2*xsize-1)*b-1)
        end
    end
    
    if in(lowercase(pbc), Set(["b", "all"]))
        for a in 1:xsize
            nfrags -= 1
            push!(chains, nfrags)
        end
    end
    
    if chains == Int64[]
        error("The PBC $pbc is not implemented yet.")
    else
        return chains
    end
end

"""
    lattice2basis(lattice::Vector{Int64}, parameters::Vector{Vector{Float64}})

Convert the lattice vector to the basis vectors of the unit cell on the cartesian space.
"""
function lattice2basis(lattice::Vector{Int64}, parameters::Vector{Vector{Float64}})
    
    a, b, c = parameters[1][1], parameters[1][2], parameters[1][3]
    α, β, γ = parameters[2][1], parameters[2][2], parameters[2][3]
    
    basis1 = [
        lattice[1]*a,
        0.,
        0.
    ]
    
    basis2 = [
        lattice[2]*b*cosd(γ),
        lattice[2]*b*sind(γ),
        0.
    ]
    
    basis3 = [
        lattice[3]*c*cosd(β),
        lattice[3]*c*(cosd(α) - cosd(β)*cosd(γ))/sind(γ),
        lattice[3]*c*sqrt(1 - cosd(β)^2 - ((cosd(α)-cosd(β)*cosd(γ))/sind(γ))^2)
    ]
    
    return [ basis1, basis2, basis3 ]
end

"""
    lattice2basis(lattice::Vector{Int64}, phase::String)
"""
function lattice2basis(lattice::Vector{Int64}, phase::String)
    return lattice2basis(lattice, get_crystallographic_info(phase)[3])
end

function firstpick(label::String)
    if isempty(label)
        return ""
    else
        return string(label[1])
    end
end

function picking_fragments(vmdoutput::Vector{SubString{String}})
    for line in vmdoutput
        if occursin(" Fragments: ", line)
            return parse(Int64, split(line)[3])
        end
    end
    error("The VMD output file does not contain the number of fragments.")
end

function picking_fragments(xyzname::String; vmd="vmd")
    tcl = tempname() * ".tcl"
    open(tcl, "w") do file
        println(file, """
        mol new "$xyzname"
        exit
        """)
    end 
    return picking_fragments(split(execVMD(vmd, tcl), "\n"))
end

"""
    isinverted(pdbname::String; vmd="vmd")

Check if the sugar monomer is inverted or not on two atoms. The default return checks if the first monomer of the cellulose chain is upside down.
"""
function isinverted(
                pdbname::String;
                resid::Int64=1, first::String="C5", last::String="O5", axis::String="z",
                vmd="vmd", vmdDebug=false
            )
    
    tcl = tempname() * ".tcl"

    open(tcl, "w") do file
        println(file, """
        mol new $pdbname
        
        set sel1 [[atomselect top "resid $resid and name $first"] get $axis]
        set sel2 [[atomselect top "resid $resid and name $last"] get $axis]

        puts "COORDS: \$sel1 \$sel2"
        
        exit
        """)
    end

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $tcl`, String), "\n")

    if vmdDebug
        return vmdoutput
    end

    for line in vmdoutput
        if occursin("COORDS:", line)
            pos = parse.(Float64, split(line)[2:3])
            return pos[2] < pos[1] ? true : false
        end
    end

    error("The VMD output file does not contain the inverting code.")
end

function printching!(file::IOStream; id=1, nresids=1, invert=false, phase="Iβ", covalent=true)
    
    inversion = invert && lowercase(phase) == "ii"

    segid, resid = "M$id", nresids
    while resid > 1
        presid = resid - 1
        inversion ? println(file, "patch 14bb $segid:$presid $segid:$resid") : println(file, "patch 14bb $segid:$resid $segid:$presid")
        
        resid -= 1
        if resid == 1 && covalent
            inversion ? println(file, "patch 14bb $segid:$nresids $segid:$resid") : println(file, "patch 14bb $segid:$resid $segid:$nresids")
        end
    end
end

function cleaning_tmpfiles(filename::String; destination_path=pwd(), new_filename="cellulose")
    
    mv(filename * ".xyz", joinpath(destination_path, "$new_filename.xyz"), force=true)
    mv(filename * ".psf", joinpath(destination_path, "$new_filename.psf"), force=true)
    mv(filename * ".pdb", joinpath(destination_path, "$new_filename.pdb"), force=true)
    mv(filename * ".tcl", joinpath(destination_path, "$new_filename.tcl"), force=true)

    files = readdir(tempdir())
    for file in files
        if occursin(".pdb", file) || occursin(".xyz", file) || occursin(".tcl", file)
            rm(joinpath(tempdir(), file))
        end
    end

    return nothing
end

"""
    hetatm!(pdbname::String)

Change the ATOM to HETATM on the PDB file.
"""
function hetatm!(pdbname::String)
    lines = readlines(pdbname)
    for (i, line) in enumerate(lines)
        if occursin("ATOM", line)
            lines[i] = replace(line, "ATOM  " => "HETATM")
        end
    end
    open(pdbname, "w") do file
        for line in lines
            println(file, line)
        end
    end
end