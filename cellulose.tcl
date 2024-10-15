package require psfgen 
topology /tmp/jl_1cBw9CZsFS.rtf 
 
segment M1 { 
    pdb /tmp/jl_R8JKNaFKVg_1.pdb 
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
patch 14bb M1:1 M1:10 
 
regenerate angles dihedrals 
 
coordpdb /tmp/jl_R8JKNaFKVg_1.pdb M1 
 
guesscoord 
 
writepsf /tmp/cellulose.psf 
writepdb /tmp/cellulose.pdb 
