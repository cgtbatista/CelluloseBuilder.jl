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
    elseif lowercase(phase) in "ii"
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