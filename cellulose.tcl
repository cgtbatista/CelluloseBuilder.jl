package require psfgen 
topology /tmp/jl_CZk7yJ8tuR.rtf 
 
segment M1 { 
    pdb /tmp/jl_cFGzCdQ9OJ_1.pdb 
} 
patch 14bb M1:8 M1:7 
patch 14bb M1:7 M1:6 
patch 14bb M1:6 M1:5 
patch 14bb M1:5 M1:4 
patch 14bb M1:4 M1:3 
patch 14bb M1:3 M1:2 
patch 14bb M1:2 M1:1 
patch 14bb M1:1 M1:8 
segment M2 { 
    pdb /tmp/jl_cFGzCdQ9OJ_4.pdb 
} 
patch 14bb M2:8 M2:7 
patch 14bb M2:7 M2:6 
patch 14bb M2:6 M2:5 
patch 14bb M2:5 M2:4 
patch 14bb M2:4 M2:3 
patch 14bb M2:3 M2:2 
patch 14bb M2:2 M2:1 
patch 14bb M2:1 M2:8 
 
regenerate angles dihedrals 
 
coordpdb /tmp/jl_cFGzCdQ9OJ_1.pdb M1 
coordpdb /tmp/jl_cFGzCdQ9OJ_4.pdb M2 
 
guesscoord 
 
writepsf /tmp/cellulose.psf 
writepdb /tmp/cellulose.pdb 
