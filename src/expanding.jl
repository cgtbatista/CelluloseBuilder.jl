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
