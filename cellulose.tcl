package require psfgen 
topology /tmp/jl_StOTSB41c9.rtf 
 
segment M1 { 
    pdb /tmp/jl_2BSzOfY3sX_1.pdb 
} 
patch 14bb M1:10 M1:9 
patch 14bb M1:9 M1:8 
patch 14bb M1:8 M1:7 
patch 14bb M1:7 M1:6 
patch 14bb M1:6 M1:5 
patch 14bb M1:5 M1:4 
patch 14bb M1:4 M1:3 
patch 14bb M1:3 M1:2 
patch 14bb M1:2 M1:1 
 
regenerate angles dihedrals 
 
coordpdb /tmp/jl_2BSzOfY3sX_1.pdb M1 
 
guesscoord 
 
writepsf /tmp/cellulose.psf 
writepdb /tmp/cellulose.pdb 
