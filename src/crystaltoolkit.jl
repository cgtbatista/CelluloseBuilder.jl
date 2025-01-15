struct XYZ
    atoms::Vector{String}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

struct XYZs
    atoms::Vector{Vector{String}}
    x::Vector{Vector{Float64}}
    y::Vector{Vector{Float64}}
    z::Vector{Vector{Float64}}
end

struct UnitCell
    atoms::Vector{String}                                   ##
    coords::Vector{Vector{Float64}}                         ## fractional coordinates [ [x0, y0, z0], [x1, y1, z1], [x2, y2, z2], ..., [xN, yN, zN] ]
    parameters::Tuple{Vector{Float64}, Vector{Float64}}     ## [ [a, b, c], [α, β, γ] ]
end


"""
    _expanding_z(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, zsize::Int64; phase="Iβ")

This function is able to propagate the crystalline system across the z-axis.
"""
function _expanding_z(xyz::XYZ, n::Int64; phase="Iβ")

    atoms, x, y, z = String[], Float64[], Float64[], Float64[]

    if !(lowercase(phase) in Set(["ia", "iα", "ib", "iβ", "ii", "iii", "iii_i", "iiii"]))
        error("The phase $phase is not implemented yet.")
    end

    if !(lowercase(phase) in Set(["ia", "iα"]))        
        parameters = get_crystallographic_info(phase)[3]
        c = parameters[1][3]
        if lowercase(phase) in ["ib", "iβ"]
            for k in collect(1:1:n)
                append!(atoms, xyz.atoms); append!(x, xyz.x); append!(y, xyz.y); append!(z, xyz.z .+ c*(k-1));
            end
        end
        if lowercase(phase) == "ii"
            for k in collect(2:1:max_k)
                append!(atoms, xyz.atoms); append!(x, xyz.x); append!(y, xyz.y); append!(z, xyz.z .+ c*(k-1));
            end
        end
        if lowercase(phase) in Set(["iii", "iii_i", "iiii"])
            for k in collect(1:1:n)
                append!(atoms, xyz.atoms); append!(x, xyz.x); append!(y, xyz.y); append!(z, xyz.z .+ c*(k-1));
            end
        end
    else
        atoms, x, y, z = xyz.atoms, xyz.x, xyz.y, xyz.z
    end

    return XYZ(atoms, x, y, z)
end

function _trimming_xy(xyz::XYZs, dim::Vector{Int64}; phase="Iβ")

    atoms, x, y, z = String[], Float64[], Float64[], Float64[]

    if !(lowercase(phase) in Set(["ia", "iα", "ib", "iβ", "ii", "iii", "iii_i", "iiii"]))
        error("The phase $phase is not implemented yet.")
    else
        xsize, ysize = dim[1], dim[2]
        units = 1
    end

    if lowercase(phase) in Set(["ib", "iβ"])
        for j in collect(1:1:ysize), i in collect(1:1:xsize)
            at, xtemp, ytemp, ztemp = xyz.atoms[units], xyz.x[units], xyz.y[units], xyz.z[units]
            _atomselect_indexes = (eachindex(at) .== 0)
            if ((i == 1) && (j != ysize)) || ((i != 1) && (i != xsize) && (j == 1))
                _atomselect_indexes = (eachindex(at) .>= 64) .& (eachindex(at) .<= 84)
            end
            if ((j != 1) && (i == xsize)) || ((i != 1) && (i != xsize) && (j == ysize))
                _atomselect_indexes = (eachindex(at) .>= 22) .& (eachindex(at) .<= 42)
            end
            if ((i == 1) && (j == ysize)) || ((j == 1) && (i == xsize))
                _atomselect_indexes = ((eachindex(at) .>= 22) .& (eachindex(at) .<= 42)) .| ((eachindex(at) .>= 64) .& (eachindex(at) .<= 84))
            end

            units += 1
            append!(atoms, at[.!(_atomselect_indexes)])
            append!(x, xtemp[.!(_atomselect_indexes)])
            append!(y, ytemp[.!(_atomselect_indexes)])
            append!(z, ztemp[.!(_atomselect_indexes)])
        end
    end
    
    if lowercase(phase) == "ii"
        for j in collect(1:1:ysize), i in collect(1:1:xsize)
            at, xtemp, ytemp, ztemp = xyz.atoms[units], xyz.x[units], xyz.y[units], xyz.z[units]
            _atomselect_indexes = (eachindex(at) .== 0)
            if ((i == 1) && (j != ysize)) || ((i != 1) && (i != xsize) && (j == 1))
                _atomselect_indexes = (eachindex(at) .>= 55) .& (eachindex(at) .<= 72)
            end
            if ((j != 1) && (i == xsize)) || ((i != 1) && (i != xsize) && (j == ysize))
                _atomselect_indexes = (eachindex(at) .>= 19) .& (eachindex(at) .<= 36)
            end
            if ((i == 1) && (j == ysize)) || ((j == 1) && (i == xsize))
                _atomselect_indexes = ((eachindex(at) .>= 19) .& (eachindex(at) .<= 36)) .| ((eachindex(at) .>= 55) .& (eachindex(at) .<= 72))
            end
            
            units += 1
            append!(atoms, at[.!(_atomselect_indexes)])
            append!(x, xtemp[.!(_atomselect_indexes)])
            append!(y, ytemp[.!(_atomselect_indexes)])
            append!(z, ztemp[.!(_atomselect_indexes)])
        end        
    end

    if lowercase(phase) in Set(["ia", "iα"])
        for (at, xtemp, ytemp, ztemp) in zip(xyz.atoms, xyz.x, xyz.y, xyz.z)
            append!(atoms, at); append!(x, xtemp); append!(y, ytemp); append!(z, ztemp);
        end
    end

    if lowercase(phase) in Set(["iii", "iii_i", "iiii"])
        for (at, xtemp, ytemp, ztemp) in zip(xyz.atoms, xyz.x, xyz.y, xyz.z)
            append!(atoms, at); append!(x, xtemp); append!(y, ytemp); append!(z, ztemp);
        end
    end
    
    return XYZ(atoms, x, y, z)
end

"""

    get_crystallographic_info(phase::String)

This function aims to return the crystallographic information for each cellulose polymorph. The information is related to the CHARMM atomnames, the unit cell
parameters and the fractional coordinates of the asymetric unit.

    - `atomnames::Vector{String}` is the default CHARMM names for each atom of the β-Glc residue.

    - `asymmetric_unit::Vector{Vector{Float64}}` store the fractional coordinates for the asymetric unit with the
       format [ [x0, y0, z0], [x1, y1, z1], [x2, y2, z2], ..., [xN, yN, zN] ].

    - `parameters_unit::Vector{Vector{Float64}}` have crystalographic parameters on the first vector (`a`, `b`, `c`)
       and angles on the second vector (`α`, `β`, `γ`)

## Arguments

- `phase::String`: The cellulose phase. It could be `Iβ`, `Iα`, `II` or `III`.

### Examples

```jldoctest

julia > get_crystallographic_info("III_I")

```

"""

function get_crystallographic_info(phase::String)
    ## [ [a, b, c], [α, β, γ] ]
    if phase == "Ib" || phase == "Iβ"
        atomnames = [ "C1", "H1", "C2", "H2", "C3", "H3", "C4", "H4", "C5", "H5", "C6", "H61", "H62", "O2", "O3", "O4", "O5", "O6", "HO2", "HO3", "HO6" ]
        asymmetric_unit =  [
            [ 0.014,  -0.042,   0.0433], [ 0.1382, -0.0188,  0.0598], [-0.026, -0.184,  -0.0516], [-0.1508, -0.2177, -0.0546], [ 0.040,  -0.137,  -0.1848],
            [ 0.1667, -0.1323, -0.1840], [-0.007,   0.030,  -0.2250], [-0.1316, 0.0176, -0.2425], [ 0.026,   0.159,  -0.1232], [ 0.1507,  0.1821, -0.1090],
            [-0.047,   0.319,  -0.153 ], [-0.1672,  0.2959, -0.1778], [-0.0407, 0.3872, -0.0765], [ 0.061,  -0.314,  -0.0003], [-0.027,  -0.261,  -0.2759],
            [ 0.079,   0.086,  -0.3424], [-0.056,   0.099,  -0.0053], [ 0.048,  0.403,  -0.254 ], [ 0.0366, -0.3194,  0.0920], [ 0.0475, -0.2400, -0.3517],
            [-0.0144,  0.4924, -0.2850], [ 0.533,   0.455,   0.3043], [ 0.6577, 0.4616,  0.3200], [ 0.475,   0.3178,  0.2094], [ 0.3480,  0.3058,  0.2052],
            [ 0.545,   0.363,   0.0765], [ 0.6716,  0.3753,  0.0793], [ 0.482,  0.526,   0.0375], [ 0.3561,  0.5139,  0.0289], [ 0.541,   0.6574,  0.1388],
            [ 0.6670,  0.6830,  0.1370], [ 0.452,   0.815,   0.1125], [ 0.3341, 0.7843,  0.0847], [ 0.4491,  0.8778,  0.1916], [ 0.528,   0.166,   0.2520],
            [ 0.486,   0.233,  -0.0081], [ 0.563,   0.586,  -0.0806], [ 0.485,  0.607,   0.2639], [ 0.542,   0.914,   0.017 ], [ 0.5035,  0.1593,  0.3448],
            [ 0.5104,  0.2733, -0.0961], [ 0.5184,  1.0275,  0.0326]
        ]
        parameters_unit = [ [ 7.784, 8.201, 10.380 ], [ 90.0, 90.0, 96.5 ] ]
    elseif phase == "Ia" || phase == "Iα"
        atomnames = [ "C1", "H1", "C2", "H2", "C3", "H3", "C4", "H4", "C5", "H5", "C6", "H61", "H62", "O4", "O2", "O3", "O5", "O6", "HO2", "HO3", "HO6" ]
        asymmetric_unit = [
            [ 0.254,  -0.054,   0.031 ], [0.1973, -0.1585, -0.1140], [ 0.193,  -0.143,   0.234 ], [ 0.2550, -0.0383,  0.3792], [ 0.022, -0.174,   0.114],
            [-0.0404, -0.2902, -0.0201], [0.000,   0.035,  -0.003 ], [ 0.0464,  0.1410,  0.1359], [ 0.079,   0.138,  -0.175 ], [ 0.0220, 0.0485, -0.3333],
            [ 0.092,   0.374,  -0.236 ], [0.1302,  0.4540, -0.0814], [ 0.1664,  0.4464, -0.3072], [-0.163,   0.003,  -0.152 ], [ 0.211, -0.346,   0.313],
            [-0.031,  -0.243,   0.307 ], [0.239,   0.152,  -0.044 ], [-0.059,   0.371,  -0.413 ], [ 0.326,  -0.296,   0.40  ], [-0.117, -0.188,   0.235],
            [-0.039,   0.43,   -0.559 ], [0.766,   0.048,  -0.026 ], [ 0.8462,  0.1471,  0.1231], [ 0.650,   0.147,  -0.218 ], [ 0.5732, 0.0483, -0.3680],
            [ 0.566,   0.180,  -0.086 ], [0.6442,  0.2833,  0.0593], [ 0.483,  -0.038,   0.008 ], [ 0.3980, -0.1353, -0.1395], [ 0.596, -0.148,   0.175],
            [ 0.6716, -0.0623,  0.3350], [0.512,  -0.382,   0.231 ], [ 0.4236, -0.4600,  0.0751], [ 0.5849, -0.4562,  0.2993], [ 0.416, -0.009,   0.157],
            [ 0.740,   0.352,  -0.288 ], [0.455,   0.270,  -0.261 ], [ 0.680,  -0.158,   0.046 ], [ 0.457,  -0.387,   0.409 ], [ 0.839,  0.333,  -0.23],
            [ 0.422,   0.302,  -0.148 ], [0.481,  -0.508,   0.51  ]
        ]
        parameters_unit = [ [ 10.400, 6.717, 5.962 ], [ 80.37, 118.08, 114.80 ] ]
    elseif phase == "II"
        atomnames = [ "C1", "H1", "C2", "H2", "C3", "H3", "C4", "H4", "C5", "H5", "C6", "H61", "H62", "O4", "O2", "O3", "O5", "O6" ]
        asymmetric_unit = [
            [-0.043,   0.007,   0.381 ], [-0.1257, -0.1115,  0.3924], [-0.125,   0.086,   0.286 ], [-0.0411,  0.2049,  0.2764], [-0.151,  -0.003,   0.156 ],
            [-0.2396, -0.1205,  0.1677], [ 0.034,   0.008,   0.112 ], [ 0.1198,  0.1243,  0.0929], [ 0.118,  -0.057,   0.216 ], [ 0.0297, -0.1731,  0.2320],
            [ 0.298,  -0.053,   0.183 ], [ 0.2935, -0.0928,  0.0947], [ 0.3960,  0.0603,  0.1874], [ 0.011,  -0.091,  -0.001 ], [-0.299,   0.062,   0.334 ], 
            [-0.224,   0.069,   0.066 ], [ 0.133,   0.034,   0.333 ], [ 0.337,  -0.155,   0.270 ], [ 0.468,   0.522,  -0.150 ], [ 0.5652,  0.6374, -0.1577],
            [ 0.318,   0.506,  -0.054 ], [ 0.2243,  0.3886, -0.0506], [ 0.398,   0.559,   0.082 ], [ 0.4800,  0.6792,  0.0821], [ 0.509,   0.468,   0.119 ],
            [ 0.4232,  0.3516,  0.1370], [ 0.640,   0.474,   0.011 ], [ 0.7246,  0.5920, -0.0036], [ 0.757,   0.392,   0.037 ], [ 0.8256,  0.4356,  0.1166],
            [ 0.6794,  0.2739,  0.0479], [ 0.621,   0.539,   0.232 ], [ 0.230,   0.603 , -0.092 ], [ 0.250,   0.523,   0.170 ], [ 0.541,   0.413,  -0.108 ],
            [ 0.886,   0.418,  -0.067 ]  
        ]
        parameters_unit = [ [ 8.10, 9.03, 10.31 ], [ 90.0, 90.0, 117.1 ] ]
    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"
        atomnames = [ "C1", "C2", "C3", "C4", "C5", "C6", "O2", "O3", "O4", "O5", "O6", "H1", "H2", "H3", "H4", "H5", "H61", "H62", "HO2", "HO3", "HO6" ]
        asymmetric_unit = [
           [ 0.050961,  0.057262, 0.380527], [ 0.191389,  0.196135, 0.287129], [ 0.047193,  0.160362, 0.153425], [ 0.009263, -0.031267,  0.112801],
           [-0.148745, -0.157932, 0.218384], [-0.176129, -0.348168, 0.188882], [ 0.181636,  0.366971, 0.331722], [ 0.242720,  0.281823,  0.066537],
           [-0.191038, -0.076143, 0.000711], [ 0.038639, -0.114787, 0.333349], [-0.238242, -0.451795, 0.303957], [-0.164524,  0.062564,  0.393673],
           [ 0.411219,  0.197152, 0.278760], [-0.158495,  0.183912, 0.154995], [ 0.213656, -0.050262, 0.093348], [-0.356184, -0.141547,  0.235149],
           [-0.342927, -0.390800, 0.126975], [ 0.016046, -0.359921, 0.150091], [ 0.40,      0.442,    0.33    ], [ 0.15,      0.27,     -0.020   ],
           [-0.12,     -0.54,     0.30]
        ]
        parameters_unit = [ [ 4.450, 7.850, 10.310 ], [ 90.0, 90.0, 105.10 ] ]
    else
        error("I'm Sorry, but the phase $phase is not available...")
    end
    return atomnames, asymmetric_unit, parameters_unit
end