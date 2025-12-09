"""
    fractional2cartesian(unitcell::Vector{Int64}, phase::String)

Convert the fractional unit cell coordinates to the cartesian coordinates. The return is the `x`, `y`, and `z` cartesian coordinates.
"""
function fractional2cartesian(unitcell::Vector{T}, phase::String) where T
    phase = lowercase(phase)
    coords = translate(phase)
    parameters = get_crystallographic_info(phase)[3]
    alpha = in(phase, ["ia", "iα"]) ? true : false
    return fractional2cartesian(
        unitcell, coords, parameters, istriclinic=alpha
    )
end

function fractional2cartesian(
    unitcell::Vector{Int64},
    fcoords::Vector{Vector{Float64}},
    parameters::Vector{Vector{Float64}}; istriclinic=true
)        
    x, y, z = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    a, b, c = parameters[1][1], parameters[1][2], parameters[1][3]
    cosα, cosβ, cosγ = cosd(parameters[2][1]), cosd(parameters[2][2]), cosd(parameters[2][3])
    sinγ = sind(parameters[2][3])
    V = a * b * c * sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)

    function cartesian(xfrac, yfrac, zfrac)
        xcart = a * xfrac + b * yfrac * cosγ + c * zfrac * cosβ
        ycart = b * yfrac * sinγ + c * zfrac * (cosα - cosβ * cosγ) / sinγ
        zcart = V * zfrac / (a * b * sinγ)
        return xcart, ycart, zcart
    end

    imax, jmax, kmax = if istriclinic
            unitcell[3]-1, unitcell[2]-1, unitcell[1]-1 ## parallelepiped
        else
            unitcell[1]-1, unitcell[2]-1, 0
    end
    
    for k in 0:kmax, j in 0:jmax, i in 0:imax
        xtemp, ytemp, ztemp = Float64[], Float64[], Float64[]
        for fcoord in fcoords
            xcart, ycart, zcart = cartesian(fcoord[1]+i, fcoord[2]+j, fcoord[3]+k)
            push!(xtemp, xcart); push!(ytemp, ycart); push!(ztemp, zcart)
        end
        push!(x, xtemp); push!(y, ytemp); push!(z, ztemp)
    end

    return x, y, z
end

function atomnames(labels::Vector{String}, n::Int64, ncells::Int64)
    atoms = repeat(replace.(labels, r"\d" => ""), outer=n);
    return [ atoms for _ in 1:ncells ]
end

function atomnames(phase::String; ncells=1, original=true)
    
    labels = get_crystallographic_info(phase)[1]

    if lowercase(phase) in ["ib", "iβ"]
        n = 4
    elseif lowercase(phase) in ["ia", "iα"]
        n = 2
    elseif lowercase(phase) == "ii"
        n = 4
    elseif lowercase(phase) in ["iii", "iii_i", "iiii"]
        n = 2
    else
        error("The phase $phase is not recognized.")
    end
    
    if original
        return atomnames(firstpick.(labels), n, ncells), labels
    else
        return atomnames(firstpick.(labels), n, ncells)
    end
end

#### Future adaptations

@inline function atomname2element(atomname::String)
    if isempty(atomname)
        throw(ArgumentError("The atom name $atomname is not valid to pass as an element."))
    end
    element = atomname[1]
    if !in(uppercase(element), Set(["C", "O", "H", "D", "N"]))
        throw(ArgumentError("The atom name $atomname is not valid to pass as an element."))
    end
    return element
end

"""
    translate(coords::Vector{Vector{Float64}}, phase::String)

Transforms the coordinates of the asymmetric unit cell to attend the translational symmetry of the cellulose phase (`Iα`, `Iβ`, `II` or `III`).
This function returns the fractional coordinates of the unit cell.
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
