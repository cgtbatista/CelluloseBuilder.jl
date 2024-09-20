"""
    get_residuePDB(residue::String; pdbname=nothing)
    
    This function creates a PDB file for a given residue.

#### Arguments
- `residue`: The residue for which the PDB file is to be created. Currently, only `PETN`.
- `pdbname`: The name of the PDB file to be created. If `nothing`, it will generate a temporary file
"""
function get_residuePDB(residue::String; pdbname=nothing)
    
    if isnothing(pdbname)
        pdbname = tempname() * ".pdb"
    end

    pdb = Base.open(filename, "w")

    Base.write(pdb, "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")

    if residue == "PETN"
        Base.write(pdb, "ATOM      1  N   PETNU   1       6.758  -1.101   4.264  1.00  0.00      U    N\n")
        Base.write(pdb, "ATOM      2  HN1 PETNU   1       7.717  -0.971   3.936  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      3  HN2 PETNU   1       6.637  -0.657   5.194  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      4  HN3 PETNU   1       6.023  -0.537   3.647  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      5  C2  PETNU   1       6.341  -2.548   4.259  1.00  0.00      U    C\n")
        Base.write(pdb, "ATOM      6  H2A PETNU   1       7.200  -3.148   3.946  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      7  H2B PETNU   1       6.059  -2.824   5.275  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM      8  C1  PETNU   1       5.171  -2.812   3.308  1.00  0.00      U    C\n")
        Base.write(pdb, "ATOM      9  H1A PETNU   1       5.096  -3.895   3.175  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM     10  H1B PETNU   1       5.368  -2.353   2.333  1.00  0.00      U    H\n")
        Base.write(pdb, "ATOM     11  P   PETNU   1       3.704  -0.781   4.106  1.00  0.00      U    P\n")
        Base.write(pdb, "ATOM     12  O3  PETNU   1       2.278  -0.388   4.235  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     13  O4  PETNU   1       4.698  -0.069   3.168  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     14  O1  PETNU   1       3.959  -0.457   5.017  0.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     15  O2  PETNU   1       3.901  -2.400   3.806  1.00  0.00      U    O\n")
        Base.write(pdb, "ATOM     16  HO  PETNU   1       4.785  -0.937   5.314  0.00  0.00      U    H\n")
    else
        throw(ArgumentError("The residue $residue is still not covered by this script."))
    end

    Base.write(pdb, "END\n")

    Base.close(pdb)

    return pdbname

end

function get_petncellulose()
    return nothing
end


