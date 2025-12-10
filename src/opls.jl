struct CelluloseFiles
    psf::AbstractString
    pdb::AbstractString
    xyz::AbstractString
    tcl::AbstractString
end

const charmm2oplsaa = Dict{Tuple{String, String}, String}([
    [("C1","CC3162")]                                                                .=> "CO" ;
    [("C2","CC3161"), ("C3","CC3161"), ("C4","CC3161")]                              .=> "CT" ;
    [("C5","CC3163")]                                                                .=> "CT" ;
    [("C6","CC321")]                                                                 .=> "CT" ;
    [("O1","OC311"), ("O2","OC311"), ("O3","OC311"), ("O4","OC311")]                 .=> "OH" ;
    [("O6","OC311")]                                                                 .=> "OH2";
    [("O5","OC3C61"), ("O4","OC301")]                                                .=> "OS" ;
    [("HO1","HCP1"), ("HO2","HCP1"), ("HO3","HCP1"), ("HO4","HCP1"), ("HO6","HCP1")] .=> "HO" ;
    [("H1","HCA1"), ("H2","HCA1"), ("H3","HCA1"), ("H4","HCA1"), ("H5","HCA1")]      .=> "HC" ;
    [("H61","HCA2"),("H62","HCA2")]                                                  .=> "HC" ;
])

function shambles_atomtype(key::Tuple{String, String})
    atomtype = get(charmm2oplsaa, key, nothing)
    if isnothing(atomtype)
        return nothing
    end
    spacing = length(key[2]) - length(atomtype)
    return rpad(atomtype, spacing)
end

function atomtypesPSF!(psfname::AbstractString)
    Base.open(psfname, "r+") do file
        lines = readlines(file)
        for (i, line) in enumerate(lines)
            if occursin("!NATOM", line)
                natoms = parse(Int, split(line)[1])
                for j in 1:natoms
                    atomline = lines[i + j]
                    name   = string(strip(split(atomline)[5]))
                    charmm = string(strip(split(atomline)[6]))
                    oplsaa = shambles_atomtype(
                        (name, charmm)
                    )
                    if !isnothing(oplsaa)
                        newline = replace(atomline, charmm => oplsaa)
                        lines[i + j] = newline
                    end
                end
                break
            end
        end
        seekstart(file)
        write(
            file, join(lines, "\n")
        )
    end
end

#
# ### ORCA related functions
#

function write_orca_input(
    pdbfile::AbstractString;
    selection="all", level="RKS RIJCOSX M062X cc-VDZ TightSCF", charge=0, multiplicity=1,
    population="hirshfeld", memory=4000, nprocs=4, filename=nothing
)
    atoms = PDBTools.readPDB(pdbfile; selection=selection)
    return write_orca_input(
        atoms;
        level=level, charge=charge, multiplicity=multiplicity,
        population=population, memory=memory, nprocs=nprocs, filename=filename
    )
end

function write_orca_input(
    atoms::Vector{<:PDBTools.Atom};
    level="RKS RIJCOSX M062X cc-pVDZ TightSCF", charge=0, multiplicity=1,
    population="hirshfeld", memory=4000, nprocs=4, filename=nothing
)
    filename = isnothing(filename) ? tempname() * ".inp" : filename
    maxcore = Int(memory / nprocs)
    open(filename, "w") do io
        println(io, """! $level PMODEL
        %output
        Print[P_$population] 1
        end
        %maxcore $maxcore
        %pal nprocs $nprocs end
        *xyz $charge $multiplicity
        """)
        for atom in atoms
            element = _opls_aa_bglc[atom.name].element
            #@printf(io, "%-2s %12.6f %12.6f %12.6f\n", element, atom.x, atom.y, atom.z)
        end
        println(io, "*")
    end
end

const orca_not_found = """
ORCA executable not found. Please, check if ORCA is installed or in your PATH/LD_LIBRARY_PATH.

You can check ORCA versions at: https://orcaforum.kofo.mpg.de/app.php/portal
"""
function execORCA(exec::String, input::String)
    hasExecutable(exec) || error(orca_not_found)
    isfile(input) || error("File not found: $input. Be sure to provide a readable file too.")
    try
        cmd = `$(exec) $input`
        return read(cmd, String)
    catch e
        throw(ErrorException("Error while executing ORCA: $e"))
    end
end

#
# ### CM5
#

function orca_hirshfeld_charges(logfile::AbstractString, natoms::Int64)
    isfile(logfile) || error("File not found: $logfile. Be sure to provide a readable file too.")
    q = Float64[]
    lines = readlines(logfile)
    for (i, line) in enumerate(lines)
        if !occursin("HIRSHFELD ANALYSIS", line)
            continue
        end
        for j in (i+7):(i+6+natoms)
            qj = parse(Float64, split(lines[j])[3])
            push!(q, qj)
        end
        return q
    end
end

const covalent_radii = Dict{String, Float64}(
    "C" => 0.760,
    "O" => 0.660,
    "H" => 0.310,
)

const D = Dict{Tuple{String, String}, Float64}(
    ("H","C") =>  0.0502,
    ("C","H") => -0.0502,
    ("H","O") =>  0.1671,
    ("O","H") => -0.1671,
    ("C","O") =>  0.0234,
    ("O","C") => -0.0234,
)

function q_cm5(atoms::Vector{<:PDBTools.Atom}, q::Vector{Float64}; α=2.474)
    q_cm5 = copy(q)
    for i in eachindex(atoms, q)
        qi = q[i]
        for j in eachindex(atoms)
            i == j && continue
            iatom, jatom = atoms[i], atoms[j]
            A, B = PDBTools.element(iatom), PDBTools.element(jatom)
            Tij = get(D, (A, B), 0.0000)                                               ## pairwise parameter
            ri = SVector{3}(iatom.x, iatom.y, iatom.z)
            rj = SVector{3}(jatom.x, jatom.y, jatom.z)
            Pij = exp(
                -α * (LinearAlgebra.norm(ri - rj) - covalent_radii[A] - covalent_radii[B])
            )                                                                          ## Pauling bond order, originaly Bij
            qi += Tij*Pij
        end
        q_cm5[i] = qi
    end
    return q_cm5
end

function adjusting_charges(
    atomtypes::Dict{Tuple{String, String}, Vector{Int64}}, charges::Vector{Float64}; λ = 1.2
)
    adjusted_charges = Dict{Tuple{String, String}, Float64}()
    charges .-= sum(charges) / length(charges)
    for (key, indices) in atomtypes
        average_charge = sum(charges[indices]) / length(indices)
        adjusted_charges[key] = round(λ * average_charge, digits=6)
    end
    return adjusted_charges
end

function atomtypes_generator(psfname::AbstractString, pdbname::AbstractString)
    atomtypes = Dict{Tuple{String, String}, Vector{Int}}()
    Base.open(psfname, "r") do file
        indices = PDBTools.index.(
            PDBTools.readPDB(pdbname)
        )
        counter = 1
        lines = readlines(file)
        for (i, line) in enumerate(lines)
            if occursin("!NATOM", line)
                natoms = parse(Int, split(line)[1])
                for j in 1:natoms
                    atomline = lines[i + j]
                    idx = parse(Int, strip(split(atomline)[1]))
                    !in(idx, indices) && continue
                    atomname = string(strip(split(atomline)[5]))
                    atomtype = string(strip(split(atomline)[6]))
                    if haskey(atomtypes, (atomname, atomtype))
                        push!(atomtypes[(atomname, atomtype)], counter)
                    else
                        atomtypes[(atomname, atomtype)] = [ counter ]
                    end
                    counter += 1
                end
            end
        end
    end
    return atomtypes
end

#
# ### Update charges
#

# maybe with String+Int64 or String+Int64+String
function chargesPSF!(psfname::AbstractString, charges::Dict{Tuple{String,String}, Float64})
    Base.open(psfname, "r+") do file
        lines = readlines(file)
        for (i, line) in enumerate(lines)
            if occursin("!NATOM", line)
                natoms = parse(Int, split(line)[1])
                for j in 1:natoms
                    atomline = lines[i + j]
                    atomname = string(strip(split(atomline)[5]))
                    atomtype = string(strip(split(atomline)[6]))
                    q_old = lpad(
                        strip(split(atomline)[7]), 9, " "
                    )
                    q_new = lpad(
                        get(charges, (atomname, atomtype), 0.0), 9, " "
                    )
                    if !isnothing(q_cm5)
                        newline = replace(atomline, q_old => q_new)
                        lines[i + j] = newline
                    end
                end
                break
            end
        end
        seekstart(file)
        write(
            file, join(lines, "\n")
        )
    end
end