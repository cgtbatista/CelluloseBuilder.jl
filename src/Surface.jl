"""

    surface(rawfile::String; surface_method="distance", selection="all", parameters="-res 0.6 -cutoff 4.0", vmd="vmd")

Generate a surface from a VMD input file using the `volmap` plugin. The surface can be generated using different methods, such as `distance`, `density`, among others.
The default method is `distance`. The selection can be defined by the user, and the default is `all`.

The parameters can be set by the user, and the default is `-res 0.6 -cutoff 4.0`. The VMD executable can be provided as an argument whether is not on default PATH.

## Arguments

- `rawfile::String`: The name of the input file without extension. The script will look for `rawfile.pdb` and `rawfile.psf`.

- `surface_method::String="distance"`: The method to be used to generate the surface using the volmap tools. The default is `distance`.
    `density`
        creates a map of the weighted atomic density at each gridpoint. This is done by replacing each atom in the selection with a normalized gaussian distribution
        of width (standard deviation) equal to its atomic radius. The gaussian distribution for each atom is then weighted using an optional -weight argument, and
        defaults to a weight of one (i.e, the number density). The various gaussians are then additively distributed on a grid.
    `interp`
        creates a map with the atomic weights interpolated onto a grid. For each atom, its weight is distributed to the 8 nearest voxels via a trilinear interpolation.
        The optional -weight argument defaults to a weight of one.
    `distance`
        creates a map for which each gridpoint contains the distance between that point and the edge of the nearest atom. In other words, each gridpoint specifies the
        maximum radius of a sphere cnetered at that point which does not intersect with the spheres of any other atoms. All atoms are treated as spheres using the
        atoms' VMD radii.
    `coulomb`, `coulombmsm`
        Creates a map of the electrostatic field of the atom selection, made by computing the non-bonded Coulomb potential from each atom in the selection (in units
        of **k_B.T/e**). The coulomb map generation is optimized to take advantage of multi-core CPUs and programmable GPUs if they are available.
    `ils`
        Creates a free energy map of the distribution of a weakly-interacting monoatomic or diatomic gas ligand throughout the system using the Implicit Ligand
        Sampling (ILS) technique. See additional information about ILS below.
    `mask`
        Creates a map which is set to 0 or 1 depending on whether they are within a specified -cutoff distance argument of any atoms in the selection. The mask map is
        typically used in combination with other maps in order to hide/mask data that is far from a region of interest.
    `occupancy`
        Each grid point is set to either 0 or 1, depending on whether it contains onbe or more atoms or not. When averaged over many frames, this will provide the
        fractional occupancy of that grid point. By default, atoms are treated as spheres using the atomic radii and a gridpoint is considered to be "occupied" if
        it lies inside that sphere. Use the -points argument to treat atoms as points (a grid point is "occupied" if its grid cube contains an atom's center). 

- `selection::String="all"`: The selection to be used to generate the surface. The default is `all`.

- `parameters::String="-res 0.6 -cutoff 4.0"`: The parameters to be used in the `volmap` plugin. The default is `-res 0.6 -cutoff 4.0`. To understand the parameters,
   see the VMD documentation https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node157.html.

- `vmd::String="vmd"`: The VMD executable. The default is `vmd`.

### Examples

```julia-repl

julia > surface("new_crystal")

```

"""
function surface(rawfile::String; surface_method="distance", selection="all", parameters="-res 0.6 -cutoff 4.0", vmd="vmd", vmdDebug=false)

    if !isfile(rawfile * ".pdb") || !isfile(rawfile * ".psf")
        error("File not found: Be sure that exist $rawfile.psf and $rawfile.pdb file.")
    end

    tcl = tempname() * ".tcl"

    open(tcl, "w") do file
        println(file, """
        mol new "$rawfile.psf"
        mol addfile "$rawfile.pdb"
        
        [atomselect top "name C1 C2 C3 C4 C5"] set radius 2.000
        [atomselect top "name C6"] set radius 2.010
        [atomselect top "name O5"] set radius 1.650
        [atomselect top "name O1 O2 O3 O4 O6"] set radius 1.765
        [atomselect top "name H1 H2 H3 H4 H5 H61 H62"] set radius 1.340
        [atomselect top "name HO1 HO2 HO3 HO4 HO6"] set radius 1.800

        volmap $surface_method [atomselect top "$selection"] $parameters -o $rawfile.dx
        
        mol delete top
        
        mol new "$rawfile.dx"
        
        axes location off
        
        mol modstyle 0 top isosurface 0.5 0 0 0 1 1
        set unitary_space {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}
        
        molinfo top set center_matrix [list \$unitary_space]
        molinfo top set rotate_matrix [list \$unitary_space]
        molinfo top set scale_matrix [list \$unitary_space]
        
        render STL $rawfile.stl true
        
        set dat_filename [open $rawfile.dat w]
        set tmp_filename [open $rawfile.stl r]
        while {[gets \$tmp_filename line] >= 0} {
          if {[lindex [split \$line] 7] == "vertex"} {
            puts \$dat_filename "[lrange [split \$line] 8 end]"
          }
        }
        close \$tmp_filename
        close \$dat_filename
        mol delete all

        exit
        """)
    end

    vmdoutput = execVMD(vmd, tcl)

    if vmdDebug
        return vmdoutput
    else
        return rawfile * ".stl", rawfile * ".dat"
    end
end