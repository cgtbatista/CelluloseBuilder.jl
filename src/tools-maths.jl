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