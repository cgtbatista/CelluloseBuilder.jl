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
    Base.write(pkminput, "tolerance 2.0\n")
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

function fibril_segid(psfname::String, pdbname::String, old_starting::String, new_starting::String; new_pdbname=nothing, vmd="vmd")
  
    new_pdbname = isnothing(new_pdbname) ? new_pdbname = tempname() * ".pdb" : new_pdbname
    new_psfname = replace(new_pdbname, ".pdb" => ".psf")

    pdb_segnames = PDBTools.segname.(PDBTools.readPDB(pdbname))
    new_segnames = replace.(pdb_segnames, old_starting => new_starting)

    tcl = tempname() * ".tcl"

    vmdinput = Base.open(tcl, "w")

    Base.write(vmdinput, "mol new $psfname\n")
    Base.write(vmdinput, "mol addfile $pdbname\n")
    Base.write(vmdinput, "\n")
    
    for ith_segids in eachindex(new_segnames)
        Base.write(vmdinput, "[atomselect top \"index $(ith_segids-1)\"] set segid \"$(new_segnames[ith_segids])\"\n")
    end
    
    Base.write(vmdinput, "\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepsf $new_psfname\n")
    Base.write(vmdinput, "[atomselect top \"all\"] writepdb $new_pdbname\n")
    Base.write(vmdinput, "exit\n")

    Base.close(vmdinput)

    vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e $(tcl)`, String), "\n")

    return new_psfname, new_pdbname, vmdoutput

end