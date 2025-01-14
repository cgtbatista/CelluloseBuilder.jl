"""
    translate(coords::Vector{Vector{Float64}}, phase::String)

Transforms the coordinates of the asymmetric unit cell to attend the translational symmetry of the cellulose phase (`Iα`, `Iβ`, `II` or `III`). This function
returns the fractional coordinates of the unit cell.
"""
function translate(coords::Vector{Vector{Float64}}, phase::String)

    xtemp, ytemp, ztemp = Float64[], Float64[], Float64[]

    valid_phases = Set(["iβ", "ib", "iα", "ia", "ii", "iii", "iii_i", "iiii"])

    if !in(lowercase(phase), valid_phases)
        error("The phase $phase is not implemented yet.")
    end

    if lowercase(phase) in ["ia", "iα"]
        for coord in coords
            push!(xtemp, coord[1])
            push!(ytemp, coord[2])
            push!(ztemp, coord[3])
        end
    else
        for coord in enumerate(repeat(coords, outer=2))
            if coord[1] <= length(coords)
                push!(xtemp, coord[2][1])
                push!(ytemp, coord[2][2])
                push!(ztemp, coord[2][3])
            else
                push!(xtemp, -1*coord[2][1])
                push!(ytemp, -1*coord[2][2])
                push!(ztemp, coord[2][3]+0.5)
            end
        end
    end

    return [ [ i, j, k ] for (i, j, k) in zip(xtemp, ytemp, ztemp) ]
end

"""
    translate(phase::String)
"""
function translate(phase::String)
    return translate(
                get_crystallographic_info(phase)[2],
                phase
            )
end

"""
    translate(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, phase::String)
"""
function translate(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, phase::String)
    
    if length(x) != length(y) || length(y) != length(z)
        error("The number of x, y, and z coordinates must be the same.")
    end

    return translate(
                [ [ i, j, k ] for (i, j, k) in zip(x, y, z) ],
                phase
            )
end