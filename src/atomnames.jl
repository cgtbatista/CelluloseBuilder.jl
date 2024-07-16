function atomsvecString(atomnames::Vector{String}, n::Int64, xsize::Int64, ysize::Int64)
    atoms = Vector{String}[]
    labels = repeat(replace.(atomnames, r"\d" => ""), outer=n);
    for j in collect(1:1:ysize), i in collect(1:1:xsize)
        push!(atoms, labels)
    end
    return atoms
end

function atomsvecString(phase::String, xsize::Int64, ysize::Int64)
    atomnames = map(_atomsvec_pick_first, get_crystallographic_info(phase)[1])
    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"; n = 4; end
    return atomsvecString(atomnames, n, xsize, ysize), get_crystallographic_info(phase)[1]
end

function atomsvecString(phase::String)
    atomnames = map(_atomsvec_pick_first, get_crystallographic_info(phase)[1])
    if phase == "I-BETA" || phase == "Ib" || phase == "Iβ"; n = 4; end
    labels = repeat(replace.(atomnames, r"\d" => ""), outer=n);
    return labels, get_crystallographic_info(phase)[1]
end

function _atomsvec_remove_numbers(string::String)
    return replace(string, r"\d" => "")
end

function _atomsvec_pick_first(string::String)
    return "$(isempty(string) ? "" : string[1])"
end