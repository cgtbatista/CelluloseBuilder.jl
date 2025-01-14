"""
    lattice2basis(lattice::Vector{Int64}, parameters::Vector{Vector{Float64}})

Convert the lattice vector to the basis vectors of the unit cell on the cartesian space.
"""
function lattice2basis(lattice::Vector{Int64}, parameters::Vector{Vector{Float64}})
    
    a, b, c = parameters[1][1], parameters[1][2], parameters[1][3]
    α, β, γ = parameters[2][1], parameters[2][2], parameters[2][3]
    
    basis1 = [
        lattice[1]*a,
        0.,
        0.
    ]
    
    basis2 = [
        lattice[2]*b*cosd(γ),
        lattice[2]*b*sind(γ),
        0.
    ]
    
    basis3 = [
        lattice[3]*c*cosd(β),
        lattice[3]*c*(cosd(α) - cosd(β)*cosd(γ))/sind(γ),
        lattice[3]*c*sqrt(1 - cosd(β)^2 - ((cosd(α)-cosd(β)*cosd(γ))/sind(γ))^2)
    ]
    
    return [ basis1, basis2, basis3 ]
end

"""
    lattice2basis(lattice::Vector{Int64}, phase::String)
"""
function lattice2basis(lattice::Vector{Int64}, phase::String)
    parameters = get_crystallographic_info(phase)[3]
    return lattice2basis(lattice, parameters)
end

function firstpick(label::String)
    if isempty(label)
        return ""
    else
        return string(label[1])
    end
end