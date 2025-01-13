"""
    _expanding_z(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, zsize::Int64; phase="Iβ")

This function is able to propagate the crystalline system across the z-axis.
"""
function _expanding_z(atoms::Vector{String}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, zsize::Int64; phase="Iβ")

    xcoords = Float64[]; ycoords = Float64[]; zcoords = Float64[]; atomnames = String[];
    parameters = get_crystallographic_info(phase)[3]
 
    if phase == "Ib" || phase == "Iβ"

        c = parameters[1][3];
        for k in collect(1:1:zsize)
            append!(atomnames, atoms); append!(xcoords, x); append!(ycoords, y); append!(zcoords, z .+ c*(k-1));
        end

    elseif phase == "Ia" || phase == "Iα"

        atomnames, xcoords, ycoords, zcoords = atoms, x, y, z

    elseif phase == "II"

        c = parameters[1][3]; max_k = zsize + 1;
        for k in collect(2:1:max_k)
            append!(atomnames, atoms); append!(xcoords, x); append!(ycoords, y); append!(zcoords, z .+ c*(k-1));
        end

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"

        c = parameters[1][3];
        for k in collect(1:1:zsize)
            append!(atomnames, atoms); append!(xcoords, x); append!(ycoords, y); append!(zcoords, z .+ c*(k-1));
        end

    else error("The phase $phase is not implemented yet."); end

    return atomnames, xcoords, ycoords, zcoords

end

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

function _trimming_xy(atoms::Vector{Vector{String}}, x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, z::Vector{Vector{Float64}}, cell_dim::Vector{Int64}; phase="Iβ")

    atomnames, xcoords, ycoords, zcoords = String[], Float64[], Float64[], Float64[]
    xsize, ysize = cell_dim[1], cell_dim[2]

    units = 1

    if phase == "Ib" || phase == "Iβ"

        for j in collect(1:1:ysize), i in collect(1:1:xsize)
            atomstemp, xtemp, ytemp, ztemp = atoms[units], x[units], y[units], z[units]
            _atomselect_indexes = (eachindex(atomstemp) .== 0)
            if ((i == 1) && (j != ysize)) || ((i != 1) && (i != xsize) && (j == 1))
                _atomselect_indexes = (eachindex(atomstemp) .>= 64) .& (eachindex(atomstemp) .<= 84)
            end
            if ((j != 1) && (i == xsize)) || ((i != 1) && (i != xsize) && (j == ysize))
                _atomselect_indexes = (eachindex(atomstemp) .>= 22) .& (eachindex(atomstemp) .<= 42)
            end
            if ((i == 1) && (j == ysize)) || ((j == 1) && (i == xsize))
                _atomselect_indexes = ((eachindex(atomstemp) .>= 22) .& (eachindex(atomstemp) .<= 42)) .| ((eachindex(atomstemp) .>= 64) .& (eachindex(atomstemp) .<= 84))
            end
            append!(atomnames, atomstemp[.!(_atomselect_indexes)])
            append!(xcoords, xtemp[.!(_atomselect_indexes)]); append!(ycoords, ytemp[.!(_atomselect_indexes)]); append!(zcoords, ztemp[.!(_atomselect_indexes)]);

            units += 1
        end
        
    elseif phase == "Ia" || phase == "Iα"

        for (at,i,j,k) in zip(atoms, x, y, z)
            append!(atomnames, at)
            append!(xcoords, i); append!(ycoords, j); append!(zcoords, k);
        end

    elseif phase == "II"

        for j in collect(1:1:ysize), i in collect(1:1:xsize)
            atomstemp, xtemp, ytemp, ztemp = atoms[units], x[units], y[units], z[units]
            _atomselect_indexes = (eachindex(atomstemp) .== 0)
            if ((i == 1) && (j != ysize)) || ((i != 1) && (i != xsize) && (j == 1))
                _atomselect_indexes = (eachindex(atomstemp) .>= 55) .& (eachindex(atomstemp) .<= 72)
            end
            if ((j != 1) && (i == xsize)) || ((i != 1) && (i != xsize) && (j == ysize))
                _atomselect_indexes = (eachindex(atomstemp) .>= 19) .& (eachindex(atomstemp) .<= 36)
            end
            if ((i == 1) && (j == ysize)) || ((j == 1) && (i == xsize))
                _atomselect_indexes = ((eachindex(atomstemp) .>= 19) .& (eachindex(atomstemp) .<= 36)) .| ((eachindex(atomstemp) .>= 55) .& (eachindex(atomstemp) .<= 72))
            end
            append!(atomnames, atomstemp[.!(_atomselect_indexes)])
            append!(xcoords, xtemp[.!(_atomselect_indexes)]); append!(ycoords, ytemp[.!(_atomselect_indexes)]); append!(zcoords, ztemp[.!(_atomselect_indexes)]);

            units += 1
        end

    elseif phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"

        for (at,i,j,k) in zip(atoms, x, y, z)
            append!(atomnames, at)
            append!(xcoords, i); append!(ycoords, j); append!(zcoords, k);
        end

    else error("The phase $phase is not implemented yet."); end
    
    return atomnames, xcoords, ycoords, zcoords

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
