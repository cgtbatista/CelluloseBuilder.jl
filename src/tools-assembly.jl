"""
   MacrofibrilAssembly(pdbname::String, fibril_radii::Float64, fibril_spacing::Float64;
   pattern="honeycomb", center=zeros(3), nlayers=nothing, outfile=nothing, old_segid_starter="M", vmd="vmd")

   Builds a macrofibril based on one elementary microfibril.
"""
function MacrofibrilAssembly(
        pdbname::String, fibril_radii::Float64, fibril_spacing::Float64;
        pattern="honeycomb", center=zeros(3), tol=2., nlayers=nothing, outfile=nothing,
        old_segid_starter="M", vmd="vmd"
    )

    if isnothing(outfile)
        outfile = tempname() * ".pdb"
    end

    nlayers = ifelse(isnothing(nlayers), 1, nlayers)  

    if pattern == "honeycomb"
        fibril_positions = honeycomb_positions(center::Vector{Float64}, fibril_radii, nlayers; spacing=fibril_spacing)
    end
    
    inp = replace(outfile, ".pdb" => ".inp")
    pkminput = Base.open(inp, "w")

    Base.write(pkminput, "## Packmol -- creating cellulose macrofibrils\n\n")
    Base.write(pkminput, "tolerance $tol\n")
    Base.write(pkminput, "filetype  pdb\n")
    Base.write(pkminput, "\n")
    Base.write(pkminput, "output    $outfile\n")

    tcl = replace(outfile, ".pdb" => ".tcl")
    vmdinput = Base.open(tcl, "w")
    
    Base.write(vmdinput, "## Psfgen -- creating a new PSF for the cellulose macrofibril\n")
    Base.write(vmdinput, "package require psfgen\n")
    Base.write(vmdinput, "\n")

    for ith_position in eachindex(fibril_positions)
        
        new_psfname, new_pdbname = fibril_segid(replace(pdbname, ".pdb" => ".psf"), pdbname, old_segid_starter, dummy_chains()[ith_position], vmd=vmd)[1:2]
        
        Base.write(pkminput, "structure $new_pdbname\n")
        Base.write(pkminput, "  number 1\n")
        Base.write(pkminput, "  center\n")
        Base.write(pkminput, "  fixed $(fibril_positions[ith_position][1]) $(fibril_positions[ith_position][2]) $(fibril_positions[ith_position][3]) 0. 0. 0.\n")
        Base.write(pkminput, "end structure\n\n")

        Base.write(vmdinput, "readpsf $new_psfname\n")

    end

    Base.write(vmdinput, "writepsf $(replace(outfile, ".pdb" => ".psf"))\n")
    Base.write(vmdinput, "exit\n")
    Base.close(vmdinput)

    Base.close(pkminput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")
    Packmol.run_packmol(inp)

    return outfile, vmdoutput
end

function SystemBoxSolvation(
        solute_pdbname::String, solvent_pdbname::String, boxlength::Float64;
        N=1, center=zeros(3), tol=2., seed=-1, outfile=nothing
    )

    if isnothing(outfile)
        outfile = tempname() * ".pdb"
    end

    inp = replace(outfile, ".pdb" => ".inp")
    pkminput = Base.open(inp, "w")

    Base.write(pkminput, "## Packmol -- creating cellulose macrofibrils\n\n")
    Base.write(pkminput, "tolerance $tol\n")
    Base.write(pkminput, "seed      $seed\n")
    Base.write(pkminput, "filetype  pdb\n")
    Base.write(pkminput, "\n")
    Base.write(pkminput, "output    $outfile\n")
    
    Base.write(pkminput, "structure $solute_pdbname\n")
    Base.write(pkminput, "  number 1\n")
    Base.write(pkminput, "  center\n")
    Base.write(pkminput, "  fixed $(center[1]) $(center[2]) $(center[3]) 0. 0. 0.\n")   ## solute
    Base.write(pkminput, "end structure\n\n")

    Base.write(pkminput, "structure $solvent_pdbname\n")
    Base.write(pkminput, "  number $N\n")
    lower_threshold = String("$(center[1] - .5 * boxlength) $(center[2] - .5 * boxlength) $(center[3] - .5 * boxlength)")
    upper_threshold = String("$(center[1] + .5 * boxlength) $(center[1] + .5 * boxlength) $(center[1] + .5 * boxlength)")
    Base.write(pkminput, "  inside box $lower_threshold, $upper_threshold\n")
    Base.write(pkminput, "end structure\n\n")

    Base.close(pkminput)

    Packmol.run_packmol(inp)

    return outfile
end


function SystemSphereSolvation(
        solute_pdbname::String, solvent_pdbname::String, radii::Float64;
        N=1, center=zeros(3), tol=2., seed=-1, outfile=nothing
    )

    if isnothing(outfile)
        outfile = tempname() * ".pdb"
    end

    inp = replace(outfile, ".pdb" => ".inp")
    pkminput = Base.open(inp, "w")

    Base.write(pkminput, "## Packmol -- creating cellulose macrofibrils\n\n")
    Base.write(pkminput, "tolerance $tol\n")
    Base.write(pkminput, "seed      $seed\n")
    Base.write(pkminput, "filetype  pdb\n")
    Base.write(pkminput, "\n")
    Base.write(pkminput, "output    $outfile\n")

    Base.write(pkminput, "structure $solute_pdbname\n")
    Base.write(pkminput, "  number 1\n")
    Base.write(pkminput, "  center\n")
    Base.write(pkminput, "  fixed $(center[1]) $(center[2]) $(center[3]) 0. 0. 0.\n")   ## solute
    Base.write(pkminput, "end structure\n\n")

    Base.write(pkminput, "structure $solvent_pdbname\n")
    Base.write(pkminput, "  number $N\n")
    sphere_threshold = String("$(center[1]) $(center[2]) $(center[3]) $(radii)")
    Base.write(pkminput, "  inside sphere $sphere_threshold\n")
    Base.write(pkminput, "end structure\n\n")

    Base.close(pkminput)

    Packmol.run_packmol(inp)

    return outfile

end


function SolvatedMacrofibril(;pdbname="initial.pdb", previous_psfname="../macrofibril.psf", topology="../../../toppar/toppar_water_ions.inp", outfile=nothing, vmd="vmd")

    if isnothing(outfile)
        outfile = tempname() * ".pdb"
    end

    tcl = replace(outfile, ".pdb" => ".tcl")

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "## VMD -- solvating cellulose macrofibrils\n\n")
    Base.write(vmdinput, "package require psfgen\n")
    
    topology = topology != Vector{String} ? [topology] : topology
    for top in topology
        Base.write(vmdinput, "topology $top\n")
    end
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "pdbalias residue HOH TIP3\n")
    Base.write(vmdinput, "readpsf  $previous_psfname\n")
    Base.write(vmdinput, "coordpdb $pdbname\n\n")

    pdbs = leftoverPDBs(pdbname)

    for pdb in pdbs
        Base.write(vmdinput, "segment $(pdb[1]) {\n")
        Base.write(vmdinput, "      auto none\n")
        Base.write(vmdinput, "      pdb  $(pdb[2])\n")
        Base.write(vmdinput, "  }\n")
        Base.write(vmdinput, "coordpdb $(pdb[2]) $(pdb[1])\n\n")
    end

    Base.write(vmdinput, "guesscoord\n\n")

    Base.write(vmdinput, "writepsf $(replace(outfile, ".pdb" => ".psf"))\n")
    Base.write(vmdinput, "writepdb $outfile\n\n")

    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return vmdoutput
end
