package require psfgen 
topology /tmp/jl_wJAbWhWbsV.rtf 
 
segment M1 { 
    pdb /tmp/jl_i1o1OjYBS7_1.pdb 
} 
patch 14bb M1:40 M1:39 
patch 14bb M1:39 M1:38 
patch 14bb M1:38 M1:37 
patch 14bb M1:37 M1:36 
patch 14bb M1:36 M1:35 
patch 14bb M1:35 M1:34 
patch 14bb M1:34 M1:33 
patch 14bb M1:33 M1:32 
patch 14bb M1:32 M1:31 
patch 14bb M1:31 M1:30 
patch 14bb M1:30 M1:29 
patch 14bb M1:29 M1:28 
patch 14bb M1:28 M1:27 
patch 14bb M1:27 M1:26 
patch 14bb M1:26 M1:25 
patch 14bb M1:25 M1:24 
patch 14bb M1:24 M1:23 
patch 14bb M1:23 M1:22 
patch 14bb M1:22 M1:21 
patch 14bb M1:21 M1:20 
patch 14bb M1:20 M1:19 
patch 14bb M1:19 M1:18 
patch 14bb M1:18 M1:17 
patch 14bb M1:17 M1:16 
patch 14bb M1:16 M1:15 
patch 14bb M1:15 M1:14 
patch 14bb M1:14 M1:13 
patch 14bb M1:13 M1:12 
patch 14bb M1:12 M1:11 
patch 14bb M1:11 M1:10 
patch 14bb M1:10 M1:9 
patch 14bb M1:9 M1:8 
patch 14bb M1:8 M1:7 
patch 14bb M1:7 M1:6 
patch 14bb M1:6 M1:5 
patch 14bb M1:5 M1:4 
patch 14bb M1:4 M1:3 
patch 14bb M1:3 M1:2 
patch 14bb M1:2 M1:1 
patch 14bb M1:1 M1:40 
segment M2 { 
    pdb /tmp/jl_i1o1OjYBS7_4.pdb 
} 
patch 14bb M2:40 M2:39 
patch 14bb M2:39 M2:38 
patch 14bb M2:38 M2:37 
patch 14bb M2:37 M2:36 
patch 14bb M2:36 M2:35 
patch 14bb M2:35 M2:34 
patch 14bb M2:34 M2:33 
patch 14bb M2:33 M2:32 
patch 14bb M2:32 M2:31 
patch 14bb M2:31 M2:30 
patch 14bb M2:30 M2:29 
patch 14bb M2:29 M2:28 
patch 14bb M2:28 M2:27 
patch 14bb M2:27 M2:26 
patch 14bb M2:26 M2:25 
patch 14bb M2:25 M2:24 
patch 14bb M2:24 M2:23 
patch 14bb M2:23 M2:22 
patch 14bb M2:22 M2:21 
patch 14bb M2:21 M2:20 
patch 14bb M2:20 M2:19 
patch 14bb M2:19 M2:18 
patch 14bb M2:18 M2:17 
patch 14bb M2:17 M2:16 
patch 14bb M2:16 M2:15 
patch 14bb M2:15 M2:14 
patch 14bb M2:14 M2:13 
patch 14bb M2:13 M2:12 
patch 14bb M2:12 M2:11 
patch 14bb M2:11 M2:10 
patch 14bb M2:10 M2:9 
patch 14bb M2:9 M2:8 
patch 14bb M2:8 M2:7 
patch 14bb M2:7 M2:6 
patch 14bb M2:6 M2:5 
patch 14bb M2:5 M2:4 
patch 14bb M2:4 M2:3 
patch 14bb M2:3 M2:2 
patch 14bb M2:2 M2:1 
patch 14bb M2:1 M2:40 
segment M3 { 
    pdb /tmp/jl_i1o1OjYBS7_7.pdb 
} 
patch 14bb M3:40 M3:39 
patch 14bb M3:39 M3:38 
patch 14bb M3:38 M3:37 
patch 14bb M3:37 M3:36 
patch 14bb M3:36 M3:35 
patch 14bb M3:35 M3:34 
patch 14bb M3:34 M3:33 
patch 14bb M3:33 M3:32 
patch 14bb M3:32 M3:31 
patch 14bb M3:31 M3:30 
patch 14bb M3:30 M3:29 
patch 14bb M3:29 M3:28 
patch 14bb M3:28 M3:27 
patch 14bb M3:27 M3:26 
patch 14bb M3:26 M3:25 
patch 14bb M3:25 M3:24 
patch 14bb M3:24 M3:23 
patch 14bb M3:23 M3:22 
patch 14bb M3:22 M3:21 
patch 14bb M3:21 M3:20 
patch 14bb M3:20 M3:19 
patch 14bb M3:19 M3:18 
patch 14bb M3:18 M3:17 
patch 14bb M3:17 M3:16 
patch 14bb M3:16 M3:15 
patch 14bb M3:15 M3:14 
patch 14bb M3:14 M3:13 
patch 14bb M3:13 M3:12 
patch 14bb M3:12 M3:11 
patch 14bb M3:11 M3:10 
patch 14bb M3:10 M3:9 
patch 14bb M3:9 M3:8 
patch 14bb M3:8 M3:7 
patch 14bb M3:7 M3:6 
patch 14bb M3:6 M3:5 
patch 14bb M3:5 M3:4 
patch 14bb M3:4 M3:3 
patch 14bb M3:3 M3:2 
patch 14bb M3:2 M3:1 
patch 14bb M3:1 M3:40 
 
regenerate angles dihedrals 
 
coordpdb /tmp/jl_i1o1OjYBS7_1.pdb M1 
coordpdb /tmp/jl_i1o1OjYBS7_4.pdb M2 
coordpdb /tmp/jl_i1o1OjYBS7_7.pdb M3 
 
guesscoord 
 
writepsf /tmp/cellulose.psf 
writepdb /tmp/cellulose.pdb 
