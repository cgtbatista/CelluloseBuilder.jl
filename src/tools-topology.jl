function generate_cellulose_topology(; filename=nothing)

    if isnothing(filename)
        filename = tempname() * ".rtf"
    end

    rtf = Base.open(filename, "w")
    
    Base.write(rtf,
        raw"""
        ! This file is the topology file needed to build cellulose crystal (CHARMM36)
     
        read rtf card append
        36 1
        AUTOGENERATE ANGLES DIHEDRALS
        !! carbon
        MASS  -1  CC3161    12.01100 C ! C2, C3, C4 CH bound to OH
        MASS  -1  CC3162    12.01100 C ! C1 (anomeric) CH bound to OH
        MASS  -1  CC3163    12.01100 C ! C5 CH bound to exocylic CH2OH
        MASS  -1  CC321     12.01100 C ! generic acyclic CH2 carbon (hexopyranose C6)
        !! hydrogen
        MASS  -1  HCA1       1.00800 H ! aliphatic proton, CH
        MASS  -1  HCA2       1.00800 H ! aliphatic proton, CH2
        MASS  -1  HCP1       1.00800 H ! polar H
        !! oxygen
        MASS  -1  OC311     15.99940 O ! hydroxyl oxygen
        MASS  -1  OC3C61    15.99940 O ! ether in six membered ring
        MASS  -1  OC301     15.99940 O ! generic linear ether
        ! DEFAults for patching FIRSt and LAST residues
        DEFA FIRS NONE LAST NONE
        AUTOGENERATE ANGLES DIHEDRALS PATCH DRUDE
        
        !! RESIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        RESI BGLC           0.000  ! 4C1 beta-D-glucose
        
        GROU                       !
        ATOM C1   CC3162    0.340  !
        ATOM H1   HCA1      0.090  !
        ATOM O1   OC311    -0.650  !
        ATOM HO1  HCP1      0.420  !
        ATOM C5   CC3163    0.110  !
        ATOM H5   HCA1      0.090  !
        ATOM O5   OC3C61   -0.400  !
        GROU                       !
        ATOM C2   CC3161    0.140  !
        ATOM H2   HCA1      0.090  !
        ATOM O2   OC311    -0.650  !
        ATOM HO2  HCP1      0.420  !
        GROU                       !
        ATOM C3   CC3161    0.140  !
        ATOM H3   HCA1      0.090  !
        ATOM O3   OC311    -0.650  !
        ATOM HO3  HCP1      0.420  !
        GROU
        ATOM C4   CC3161    0.140  !
        ATOM H4   HCA1      0.090  !
        ATOM O4   OC311    -0.650  !
        ATOM HO4  HCP1      0.420  !
        GROU
        ATOM C6   CC321     0.050  !
        ATOM H61  HCA2      0.090  !
        ATOM H62  HCA2      0.090  !
        ATOM O6   OC311    -0.650  !
        ATOM HO6  HCP1      0.420  !
        
        BOND C1   O1        C1   H1        O1   HO1       C1   O5        C1   C2
        BOND C2   H2        C2   O2        O2   HO2       C2   C3        C3   H3
        BOND C3   O3        O3   HO3       C3   C4        C4   H4        C4   O4
        BOND O4   HO4       C4   C5        C5   H5        C5   C6        C6   H61
        BOND C6   H62       C6   O6        O6   HO6       C5   O5
        !    I    J    K    L      R(IK)   T(IKJ)    PHI   T(JKL)   R(KL)
        IC   O1   C2  *C1   H1     1.3899  110.90  120.10  104.58   1.0836
        IC   O1   O5  *C1   C2     1.3899  108.62  122.10  110.88   1.5316
        IC   O2   C3  *C2   H2     1.4594  108.12 -118.78  111.06   1.1375
        IC   O2   C1  *C2   C3     1.4594  115.65 -125.60  113.28   1.4983
        IC   O3   C4  *C3   H3     1.4071  113.48  122.06  103.39   1.0895
        IC   O3   C2  *C3   C4     1.4071  108.48  124.18  109.26   1.5497
        IC   O4   C5  *C4   H4     1.3940  111.12 -110.35  108.66   1.0857
        IC   O4   C3  *C4   C5     1.3940  112.77 -129.39  115.62   1.5530
        IC   C6   O5  *C5   H5     1.5597  111.17  120.85  110.98   1.1092
        IC   C6   C4  *C5   O5     1.5597  109.90  122.92  110.30   1.4512
        IC   O6   H62 *C6   H61    1.4589  116.11 -112.93  103.57   1.1467
        IC   O6   C5  *C6   H62    1.4589  109.41 -135.95  118.22   1.0853
        IC   O5   C1   C2   C3     1.4620  110.88   57.82  113.28   1.4983
        IC   C1   C2   C3   C4     1.5316  113.28  -48.40  109.26   1.5497
        IC   C2   C3   C4   C5     1.4983  109.26   45.07  115.62   1.5530
        IC   C3   C4   C5   O5     1.5497  115.62  -49.19  110.30   1.4512
        IC   C4   C5   O5   C1     1.5530  110.30   56.36  112.12   1.4620
        IC   C5   O5   C1   C2     1.4512  112.12  -61.39  110.88   1.5316
        IC   C4   C5   C6   O6     1.5530  109.90 -177.46  109.41   1.4589
        IC   O5   C1   O1   HO1    1.4620  108.62   72.25  106.48   0.9328
        IC   C1   C2   O2   HO2    1.5316  115.65  135.41  116.81   0.9527
        IC   C2   C3   O3   HO3    1.4983  108.48  -71.46  120.86   0.9441
        IC   C3   C4   O4   HO4    1.5497  112.77   47.45  109.31   0.9911
        IC   C5   C6   O6   HO6    1.5597  109.41  -54.60  118.82   0.95210
        PATC  FIRS NONE LAST NONE
        
        !! PATCHES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! equatorial-equatorial 1->4 linkage
        ! LACTOS03, EYOCUQ01, CELLOB01
        PRES 14bb           0.02 ! (i)1->4(i-1) equatorial at C1 and equatorial at C4
        dele atom 1HO4
        dele atom 2HO1
        dele atom 2O1
        ATOM 1C4  CC3161    0.09 !
        ATOM 1O4  OC301    -0.36 !
        ATOM 2C1  CC3162    0.29 !
        BOND 1O4  2C1
        !    I    J    K    L      R(IK)   T(IKJ)    PHI   T(JKL)   R(KL)
        IC   1C3  1C4  1O4  2C1    1.5009  110.76   81.86  121.00   1.3902  ! psi
        IC   1C4  1O4  2C1  2O5    1.4560  121.00 -130.97  108.63   1.4470  ! phi
        IC   1O4  2O5 *2C1  2C2    1.3896  108.63  122.12  110.87   1.5318
        IC   2O5  1O4 *2C1  2H1    1.4470  108.63  121.92  111.32   1.0837
        
        END
        """
        )

    Base.close(rtf)

    return filename
end

function generate_petn_topology(; filename=nothing)

    if isnothing(filename)
        filename = tempname() * ".rtf"
    end

    rtf = Base.open(filename, "w")
    
    Base.write(rtf, 
        raw"""
        ! This file is the topology file needed to generate pEtN and patch it to cellulose chain (CHARMM36)
        
        read rtf card append
        36 1
        AUTOGENERATE ANGLES DIHEDRALS
        !! already exists on cgenff
        MASS  -1  HGA2       1.00800 H                          ! alphatic proton, CH2
        MASS  -1  HGP1       1.00800 H                          ! polar H
        MASS  -1  HGP2       1.00800 H                          ! polar H, +ve charge
        MASS  -1  CG321     12.01100 C                          ! aliphatic C for CH2
        MASS  -1  CG324     12.01100 C                          ! aliphatic C for CH2, adjacent to positive N (piperidine) (+)
        MASS  -1  NG3P3     14.00700 N                          ! primary NH3+, phosphatidylethanolamine
        MASS  -1  OG2P1     15.99940 O                          ! =O in phosphate or sulfate
        MASS  -1  OG303     15.99940 O                          ! phosphate/sulfate ester oxygen
        MASS  -1  OG311     15.99940 O                          ! hydroxyl oxygen
        MASS  -1  PG1       30.97380 P                          ! phosphate -1
        !! already exists on carb force field
        MASS  -1  CC321     12.01100 C                          ! aliphatic C for CH2 (e.g. CA and CB) it should be equal C6
        MASS  -1  HCA2       1.00800 H                          ! alphatic proton, CH2
        MASS  -1  HCP1       1.00800 H                          ! polar H
        MASS  -1  OC2DP     15.99940 O                          ! =O in phosphate or sulfate		    || O2L in top_all27_lipid.rtf
        MASS  -1  OC312     15.99940 O                          ! hydroxyl oxygen			            || OHL in top_all27_lipid.rtf
        MASS  -1  OC30P     15.99940 O                          ! phosphate/sulfate ester oxygen		|| OSL in top_all27_lipid.rtf
        MASS  -1  PC        30.97380 P                          ! phosphate -1
        !! new ones to phosphoethanolamine patch
        MASS  -1  CC324     12.01100 C                          ! aliphatic C for CH2 (e.g. CA and CB) it should be equal C6
        MASS  -1  HCP2       1.00800 H                          ! polar H, +ve charge
        MASS  -1  NC3P3     14.00700 N                          ! primary NH3+, phosphatidylethanolamine
        ! DEFAults for patching FIRSt and LAST residues
        DEFA FIRS NONE LAST NONE
        AUTOGENERATE ANGLES DIHEDRALS PATCH DRUDE

        !! RESIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RESI ENP            0.000 !   Phosphoethanolamine
        GROUP             ! CHARGE
        ATOM NP     NG3P3  -0.293 !
        ATOM HN1    HGP2    0.328 !
        ATOM HN2    HGP2    0.328 !
        ATOM HN3    HGP2    0.328 !
        ATOM CA     CG324   0.130 !
        ATOM HA1    HGA2    0.090 !
        ATOM HA2    HGA2    0.090 !
        ATOM CB     CG321  -0.084 !
        ATOM HB1    HGA2    0.090 !
        ATOM HB2    HGA2    0.090 !
        ATOM P      PG1     1.505 !
        ATOM O4P    OG2P1  -0.824 !
        ATOM O3P    OG2P1  -0.824 !
        ATOM O2P    OG303  -0.623 !
        ATOM O1P    OG311  -0.669 !
        ATOM HOP    HGP1    0.338 !
        BOND NP   HN1  !
        BOND NP   HN2  !
        BOND NP   HN3  !
        BOND NP   CA   !
        BOND CA   CB   !
        BOND CA   HA1  !
        BOND CA   HA2  !
        BOND CB   HB1  !
        BOND CB   HB2  !
        BOND CB   O2P  !
        BOND P    O3P  !
        BOND P    O1P  !
        BOND P    O4P  !
        BOND P    O2P  !
        BOND O1P  HOP  !
        PATC FIRS NONE LAST NONE

        !! PATCHES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        PRES PETN         -0.270       ! pEtN-BGlc linkage
        dele atom HO6                  ! residual charge = +0.420
        ATOM C6  CC321    -0.075       ! B-Glc O6 attacking the pEtN P (like a SN1 mechanism)
        ATOM O6  OC311    -0.571       ! phosphatidylethanolamine + cellulose -> pEtN-cellulose + diacylglycerol
        ATOM NP  NC3P3    -0.293 !
        ATOM HN1 HCP2      0.328 !
        ATOM HN2 HCP2      0.328 !
        ATOM HN3 HCP2      0.328 !
        ATOM CA  CC324     0.130 !
        ATOM HA1 HCA2      0.090 !
        ATOM HA2 HCA2      0.090 !
        ATOM CB  CC321    -0.077 !
        ATOM HB1 HCA2      0.090 !
        ATOM HB2 HCA2      0.090 !
        ATOM P   PC        1.501 !
        ATOM O2P OC30P    -0.569 !
        ATOM O3P OC2DP    -0.785 !
        ATOM O4P OC2DP    -0.785 !
        BOND O6 P
        BOND NP HN1  NP HN2  NP HN3  NP CA  CA CB  CA HA1  CA HA2  CB HB1  CB HB2  CB O2P  P O2P  P O3P  P O4P
        !    I     J     K     L      R(IK)  T(IKJ)   PHI    T(JKL)   R(KL) 
        IC   HN1   CA   *NP    HN2   1.0400  110.15  122.13  117.40  1.0400 ! based on the IC of the pEtN in DLiPE
        IC   HN1   CA   *NP    HN3   1.0400  110.15 -124.28  113.58  1.0400 ! toppar_all36_lipid_miscellaneous.str
        IC   HN1   NP    CA    CB    1.0400  110.15  157.44  108.52  1.4921 !
        IC   CB    NP   *CA    HA1   1.4921  108.52 -124.86  113.65  1.1110 !
        IC   HA1   NP   *CA    HA2   1.1110  113.65 -107.22  102.84  1.1110 !
        IC   NP    CA    CB    O2P   1.5241  108.52   60.11  106.54  1.4914 !
        IC   O2P   CA   *CB    HB1   1.4914  106.54 -125.10  115.49  1.1110 !
        IC   HB1   CA   *CB    HB2   1.1110  115.49 -116.55  111.62  1.1110 !
        IC   CA    CB    O2P   P     1.4921  106.54  109.74  118.42  1.5799 !
        IC   CB    O2P   P     O6    1.4914  118.42   45.68  108.09  1.5231 !
        IC   O6    O2P  *P     O3P   1.5231  108.09  119.87  104.67  1.4621 !
        IC   O6    O2P  *P     O4P   1.5231  108.09 -120.31  106.92  1.5173 !
        IC   O2P   P     O6    C6    1.5799  108.09   83.12  132.97  1.4179 !
        IC   P     O6    C6    C5    1.5231  132.97  -54.60  109.41  1.5597 !
        END
        """
    )
    
    Base.close(rtf)

    return filename
end

function generate_waterions_topology(; filename=nothing)

    if isnothing(filename)
        filename = tempname() * ".rtf"
    end

    rtf = Base.open(filename, "w")
    
    Base.write(rtf,
        raw"""
        ! This file is the topology file needed to generate water and ions (CHARMM36)

        read rtf card @app
        * Topology for water and ions
        *
        31  1

        MASS  -1  HT         1.00800 H ! TIPS3P WATER HYDROGEN
        MASS  -1  HX         1.00800 H ! hydroxide hydrogen
        MASS  -1  OT        15.99940 O ! TIPS3P WATER OXYGEN
        MASS  -1  OX        15.99940 O ! hydroxide oxygen
        MASS  -1  LIT        6.94100 LI ! Lithium ion
        MASS  -1  SOD       22.98977 NA ! Sodium Ion
        MASS  -1  MG        24.30500 MG ! Magnesium Ion
        MASS  -1  POT       39.09830 K ! Potassium Ion
        MASS  -1  CAL       40.08000 CA ! Calcium Ion
        MASS  -1  RUB       85.46780 RB ! Rubidium Ion
        MASS  -1  CES      132.90545 CS ! Cesium Ion
        MASS  -1  BAR      137.32700 BA ! Barium Ion
        MASS  -1  ZN        65.37000 ZN ! zinc (II) cation
        MASS  -1  CAD      112.41100 CD ! cadmium (II) cation
        MASS  -1  CLA       35.45000 CL ! Chloride Ion
        default first none last none
        AUTO ANGLE DIHE

        !! RESIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RESI HOH         0.000 NOANG NODIH ! tip3p water model
        GROUP
        ATOM OH2  OT     -0.834
        ATOM H1   HT      0.417
        ATOM H2   HT      0.417
        BOND OH2 H1 OH2 H2 H1 H2    ! the last bond is needed for shake
        ANGLE H1 OH2 H2             ! required
        DONOR H1 OH2
        DONOR H2 OH2
        ACCEPTOR OH2
        PATCHING FIRS NONE LAST NONE

        RESI OH       -1.00 ! hydroxide ion by adm.jr.
        GROUP
        ATOM O1 OX    -1.32
        ATOM H1 HX     0.32
        BOND O1 H1
        DONOR H1 O1
        ACCEPTOR O1

        RESI LIT       1.00 ! Lithium Ion
        GROUP
        ATOM LIT  LIT  1.00
        PATCHING FIRST NONE LAST NONE

        RESI SOD       1.00 ! Sodium Ion
        GROUP
        ATOM SOD  SOD  1.00
        PATCHING FIRST NONE LAST NONE

        RESI MG        2.00 ! Magnesium Ion
        GROUP
        ATOM MG   MG   2.00
        PATCHING FIRST NONE LAST NONE

        RESI POT       1.00 ! Potassium Ion
        GROUP
        ATOM POT   POT 1.00
        PATCHING FIRST NONE LAST NONE

        RESI CAL       2.00 ! Calcium Ion
        GROUP
        ATOM CAL  CAL  2.00
        PATCHING FIRST NONE LAST NONE

        RESI RUB       1.00 ! Rubidium Ion
        GROUP
        ATOM RUB  RUB  1.00
        PATCHING FIRST NONE LAST NONE

        RESI CES       1.00 ! Cesium Ion
        GROUP
        ATOM CES  CES  1.00
        PATCHING FIRST NONE LAST NONE

        RESI BAR       2.00 ! Barium Ion
        GROUP
        ATOM BAR  BAR  2.00
        PATCHING FIRST NONE LAST NONE

        RESI ZN2       2.00 ! Zinc (II) cation, Roland Stote
        GROUP
        ATOM ZN   ZN   2.00
        PATCHING FIRST NONE LAST NONE

        RESI CD2       2.00 ! Cadmium (II) cation
        GROUP
        ATOM CD   CAD  2.00
        PATCHING FIRST NONE LAST NONE

        RESI CLA      -1.00 ! Chloride Ion
        GROUP
        ATOM CLA  CLA -1.00
        PATCHING FIRST NONE LAST NONE

        END
        """
    )
    
    Base.close(rtf)

    return filename
end

