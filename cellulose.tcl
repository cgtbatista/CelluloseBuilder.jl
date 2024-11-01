package require psfgen 
topology /tmp/jl_rqnUvc7cpA.rtf 
 
segment M1 { 
    pdb /tmp/jl_HQJepChAtW_1.pdb 
} 
patch 14bb M1:2 M1:1 
 
regenerate angles dihedrals 
 
coordpdb /tmp/jl_HQJepChAtW_1.pdb M1 
 
guesscoord 
 
writepsf /tmp/cellulose.psf 
writepdb /tmp/cellulose.pdb 
