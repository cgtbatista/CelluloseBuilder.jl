"""
    fractional2cartesian(unitcell::Vector{Int64}, phase::String)

Convert the fractional unit cell coordinates to the cartesian coordinates. The return is the `x`, `y`, and `z` cartesian coordinates.
"""
function fractional2cartesian(unitcell::Vector{Int64}, phase::String)

    fractional_coords = transformASU(phase)
    parameters = get_crystallographic_info(phase)[3]

    invalue = if lowercase(phase) in ["ia", "iα"]
            true
        else
            false
    end

    return fractional2cartesian(unitcell, fractional_coords, parameters, inclined=invalue)
end

function fractional2cartesian(unitcell::Vector{Int64}, fractional_coords::Vector{Vector{Float64}}, parameters::Vector{Vector{Float64}}; inclined=true)
        
    x, y, z = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]

    a, b, c = parameters[1][1], parameters[1][2], parameters[1][3]
    cosα, cosβ, cosγ = cosd(parameters[2][1]), cosd(parameters[2][2]), cosd(parameters[2][3])
    sinγ = sind(parameters[2][3])
    V = a * b * c * sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)

    imax, jmax, kmax = if inclined
            unitcell[3]-1, unitcell[2]-1, unitcell[1]-1 ## parallelepiped
        else
            unitcell[1]-1, unitcell[2]-1, 0
    end

    function cartesian(xfrac, yfrac, zfrac)
        xcart = a * xfrac + b * yfrac * cosγ + c * zfrac * cosβ
        ycart = b * yfrac * sinγ + c * zfrac * (cosα - cosβ * cosγ) / sinγ
        zcart = V * zfrac / (a * b * sinγ)
        return xcart, ycart, zcart
    end
    
    for k in 0:kmax, j in 0:jmax, i in 0:imax
        xtemp, ytemp, ztemp = Float64[], Float64[], Float64[]
        for fcoord in fractional_coords
            xcart, ycart, zcart = cartesian(fcoord[1]+i, fcoord[2]+j, fcoord[3]+k)
            push!(xtemp, xcart); push!(ytemp, ycart); push!(ztemp, zcart)
        end
        push!(x, xtemp); push!(y, ytemp); push!(z, ztemp)
    end

    return x, y, z
end
