## VMD -- solvating cellulose macrofibrils

package require psfgen

topology /home/carlos/Dropbox/works/pEtNCellulose/assembly/toppar/water_ions.str

pdbalias residue HOH TIP3
readpsf  /home/carlos/Dropbox/works/pEtNCellulose/assembly/MakeSystems/chain/ccc.psf
coordpdb /home/carlos/Dropbox/works/pEtNCellulose/assembly/MakeSystems/chain/A_nacl/Accc.pdb

segment WATA {
      auto none
      pdb  /tmp/jl_B0bWHgjqxU.pdb
  }
coordpdb /tmp/jl_B0bWHgjqxU.pdb WATA

segment WATB {
      auto none
      pdb  /tmp/jl_zKjcjwV2Xe.pdb
  }
coordpdb /tmp/jl_zKjcjwV2Xe.pdb WATB

segment WATC {
      auto none
      pdb  /tmp/jl_jALQPIBubc.pdb
  }
coordpdb /tmp/jl_jALQPIBubc.pdb WATC

segment NAD {
      auto none
      pdb  /tmp/jl_LUzQ9MXOBm.pdb
  }
coordpdb /tmp/jl_LUzQ9MXOBm.pdb NAD

segment CLE {
      auto none
      pdb  /tmp/jl_N4bUbiVcT4.pdb
  }
coordpdb /tmp/jl_N4bUbiVcT4.pdb CLE

guesscoord

writepsf teste.psf
writepdb teste.pdb

exit
