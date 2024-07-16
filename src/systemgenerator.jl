function _exporting_PDBfile(phase, max_monomer, tmp, topology_file; covalent=true, vmd="vmd")

    fraglist = collect(1:1:length(tmp))
    monomer = 0; prev_monomer = 0;

    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"
            
        vmdinput = Base.open("/tmp/patching.tcl", "w")
        
        Base.write(vmdinput, "package require psfgen \n")
        Base.write(vmdinput, "topology $(topology_file) \n")
        Base.write(vmdinput, " \n")
        for frag in fraglist
            segname = "M$frag"
            Base.write(vmdinput, "segment $segname { \n")
            Base.write(vmdinput, "    pdb $(tmp[frag]) \n")
            Base.write(vmdinput, "} \n")
            monomer = max_monomer
            while monomer > 1
                prev_monomer = monomer - 1
                Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$prev_monomer \n")
                monomer -= 1
                if monomer == 1 && covalent == true
                    Base.write(vmdinput, "patch 14bb $segname:$monomer $segname:$max_monomer \n")
                end
            end
        end
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "regenerate angles dihedrals \n")
        Base.write(vmdinput, " \n")
        for frag in fraglist
            segname = "M$frag"
            Base.write(vmdinput, "coordpdb $(tmp[frag]) $segname \n")
        end
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "guesscoord \n")
        Base.write(vmdinput, " \n")
        Base.write(vmdinput, "writepsf /tmp/cellulose.psf \n")
        Base.write(vmdinput, "writepdb /tmp/cellulose.pdb \n")
        Base.close(vmdinput)
        vmdoutput = Base.split(Base.read(`$vmd -dispdev text -e /tmp/patching.tcl`, String), "\n")
    end

    return vmdoutput
end