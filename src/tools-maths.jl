"""
    R_matrix(v1, v2)

    Picks two vectors and gives the rotation matrix `R` needed to rotate the vector v2 on v1 (the mirror vector).
    If the ```||v_cross|| ≈ 0```, we do not need to rotate the residue. Otherwise, we will pick a generalize rotation matrix `R`.
"""
function R_matrix(v1, v2)

    # Normalizing the vectors
    v1 /= norm(v1)
    v2 /= norm(v2)

    # cross rotation vector
    crossvector = cross(v1, v2)
    θ = acos(
            clamp(dot(v1, v2), -1.0, 1.0)
        ) # just to ensure the angle domain


    # checking with there is need in computation of the cross matrix
    if norm(crossvector) ≈ 0
        return I(3) # return the identity
    else
        crossvector /= norm(crossvector)
    end

    # cross matrix
    CrossMatrix = [ 0 -crossvector[3] crossvector[2]; crossvector[3] 0 -crossvector[1]; -crossvector[2] crossvector[1] 0 ]
    return I(3) + sin(θ) * CrossMatrix + (1-cos(θ)) * (CrossMatrix * CrossMatrix)

end

"""
    translate_residue(coords, translation_vector)

    Translate the residue coordinates by a given directional vector.
"""
function translate_residue(coords, translation_vector)

    for xyz in eachindex(coords)
        coords[xyz] += translation_vector
    end

    return coords
end

"""
    rotate_residue(coords, reference, rotation_matrix)

    Rotate the residue coordinates by a given rotation matrix. The `ref_center` is a reference point that we will lock while rotating the other atoms.
"""
function rotate_residue(coords, ref_center, rotation_matrix)

    for xyz in eachindex(coords)
        coords[xyz] = ref_center + rotation_matrix * (coords[xyz] - ref_center)
    end

    return coords
end

"""
    honeycomb_positions(center::Vector{Float64}, r::Float64; spacing=1.)

    Generate a honeycomb lattice with 7 units around the `center` point. The `r` is the radius of the hexagon, it can be based on the fibril mean radius.
    The `spacing` flag gives the separation between them.
"""
function honeycomb_positions(center::Vector{Float64}, r::Float64; spacing=1.)

    points = []

    d_adjusted = spacing + 2 * r
    
    push!(points, [center[1], center[2], center[3]])
    push!(points, [center[1] - d_adjusted, center[2], center[3]])
    push!(points, [center[1] + d_adjusted, center[2], center[3]])
    push!(points, [center[1] - d_adjusted/2, center[2] + d_adjusted * √3/2, center[3]])
    push!(points, [center[1] + d_adjusted/2, center[2] + d_adjusted * √3/2, center[3]])
    push!(points, [center[1] - d_adjusted/2, center[2] - d_adjusted * √3/2, center[3]])
    push!(points, [center[1] + d_adjusted/2, center[2] - d_adjusted * √3/2, center[3]])

    return points

end

"""
    honeycomb_positions(center::Vector{Float64}, r::Float64, nlayers::Int64; spacing=1.)

    Generate a honeycomb lattice with `nlayers` layers around the `center` point. The `r` is the radius of the hexagon, it can be based on the fibril mean radius.
    The `spacing` flag gives the separation between them.
"""
function honeycomb_positions(center::Vector{Float64}, r::Float64, nlayers::Int64; spacing=1.)

    points = []

    d_adjusted = spacing + 2 * r

    push!(points, center)

    for n in collect(1:1:nlayers)

        for k in 0:5
            θ = k * π / 3  # angle on radii
            for m in 1:n
                # computing the coordinates
                x = center[1] + m * d_adjusted * cos(θ)
                y = center[2] + m * d_adjusted * sin(θ)
                push!(points, [x, y, center[3]])
            end
        end
    end

    return points
end


"""

    antoine_parameters(T::Float64, component::String)

    Returns the vapor pressure based on Antoine equation for a given temperature `T` and component `component` (if you want to guess the parameters).
    Commonly, the coefficients are based on a fitting using temperature is given on Celsius degree (°C).
"""
function p_antoine(T::Float64; component="water", guess=true, A=nothing, B=nothing, C=nothing)
  
    if !guess
        if isnothing(A) || isnothing(B) || isnothing(C)
            error("You need to provide all the Antoine parameters for the component.")
        end
    else
        A, B, C = antoine_parameters(T, component)
    end

    x = A - B / (C + T)
    return (10 ^ x) / 750.062
end

"""
    
        antoine_parameters(T::Float64, component::String)
    
        Returns the Antoine parameters for a given temperature `T` and component `component`.
"""
function antoine_parameters(T::Float64, component::String; initial_coefficients=[8.0, 1800.0, 240.0])

    # Antoine equation: log10(P) = A - B / (C + T)
    antoine_model(T, coefficients) = coefficients[1] .- coefficients[2] ./ (coefficients[3] .+ T)

    # Dados de temperatura (em °C) e pressão (convertida para log10(P))
    if component == "water" || component == "Water"
        if 1. <= T <= 97.98
            temperature = [
                1.00, 8.46, 15.92, 23.38, 30.84, 38.30, 45.76, 53.22, 60.68, 68.14, 75.6, 83.06, 90.52, 97.98
            ]
            pressure = [
                0.00651326, 0.0110022, 0.018010, 0.028652, 0.044401, 0.0671718, 0.099396, 0.144111, 0.20504,
                0.286684, 0.394403, 0.534498, 0.714286, 0.942169
            ] * 750.062 ## from bar to mmHg
        elseif 105.44 <= T <= 374.0
            temperature = [
                105.44, 112.90, 120.36, 127.82, 135.28, 142.74, 150.20, 157.66, 165.12, 172.58, 
                180.04, 187.50, 194.96, 202.42, 209.88, 217.34, 224.80, 232.26, 239.72, 247.18, 
                254.64, 262.10, 269.56, 277.02, 284.48, 291.94, 299.40, 306.86, 314.32, 321.78, 
                329.24, 336.70, 344.16, 351.62, 359.08, 366.54, 374.00
            ]
            pressure = [
                1.2299, 1.57724, 2.00219, 2.51746, 3.13699, 3.87597, 4.7509, 5.77952, 6.98083, 8.37506, 
                9.98367, 11.8292, 13.9355, 16.3272, 19.0302, 22.0711, 25.4776, 29.2781, 33.5016, 38.178, 
                43.3375, 49.011, 55.2296, 62.0249, 69.4288, 77.473, 86.1897, 95.6108, 105.768, 116.694, 
                128.42, 140.977, 154.397, 168.709, 183.945, 200.133, 217.304
            ] * 750.062
        else
            error("Temperature out of range.")
        end
    else
        error("Component not found.")
    end

    # Convert pressures to log10(P)
    log_pressure = log10.(pressure)

    # Ajuste dos dados usando o modelo de Antoine
    fit = LsqFit.curve_fit(antoine_model, temperature, log_pressure, initial_coefficients)

    # Parâmetros ajustados
    parameters = fit.param

    return parameters[1], parameters[2], parameters[3]
end

function V_steam(T::Float64, N::Int64; Z=.935,
        component="water", guess=true, A=nothing, B=nothing, C=nothing                  # p_antoine flags
    )

    p = p_antoine(T, component=component, guess=guess, A=A, B=B, C=C)

    N_Avogadro = 6.02214076E+23          ## mol^(-1)
    R = 8.31446261815324E-2               ## bar × L / (K · mol)

    n = Float64(N) / N_Avogadro

    V = Z * n * R * (T + 273.15) / p

    ## converting from L to Å³
    return V * 1e27
end


function N_box(V::Float64, density::Float64; MM_solvent = 18.015, solute=false, V_solute=10.)

    N_avogadro = 6.02214076E+23
    atomic_density = (density / MM_solvent) * 10^(-24) * N_avogadro ## n molecules /A³

    if !solute
        N = V * atomic_density
    else
        N = (V - V_solute) * atomic_density
    end

    return Int64(round(N))

end

function fibril_dimensions(pdbname::String, habit::String)
    
    pdb = PDBTools.readPDB(pdbname)

    if habit == "234432"
        atom1 = PDBTools.select(pdb, atom ->
                atom.segname == "M17" && atom.resnum == 1 && atom.name == "HO6")
        atom2 = PDBTools.select(pdb, atom ->
                atom.segname == "M2" && atom.resnum == 1 && atom.name == "HO2")
        min_resid = minimum(PDBTools.resnum.(pdb))
        max_resid = maximum(PDBTools.resnum.(pdb))
        atom3 = PDBTools.select(pdb, atom ->
                atom.segname == "M1" && atom.resnum == min_resid && atom.name == "HO4")
        atom4 = PDBTools.select(pdb, atom ->
            atom.segname == "M1" && atom.resnum == max_resid && atom.name == "HO1")
    else
        throw(ArgumentError("There is no implementation to the $habit fibril..."))
    end
    
    r_fibril = 0.5 * norm(PDBTools.coor(atom1)-PDBTools.coor(atom2))
    l_fibril = norm(PDBTools.coor(atom3)-PDBTools.coor(atom4))

    return round(r_fibril, digits=1), round(l_fibril, digits=1)

end