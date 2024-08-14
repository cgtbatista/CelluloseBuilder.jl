"""

    atomsvecString(atomnames::Vector{String}, n::Int64, xsize::Int64, ysize::Int64)
    atomsvecString(phase::String, xsize::Int64, ysize::Int64)
    atomsvecString(phase::String)

This function generates a vector of atomic labels based on standard CHARMM atomnames with only the atom
name. In some occasions, it can return the original atomtypes too. Carefully, because it just not work
properly on atoms with more than 1 character (e.g. Fe, Cl, Rn, etc).

### Arguments

- `atomnames::String`: vector labels with the CHARMM atomnames.
- `n::Int64`: repetition units of `atomnames` vector.
- `xsize::Int64`: the number of unit cells along x axis (`a`)
- `ysize::Int64`: the number of unit cells along y axis (`b`)
- `phase::String`: the cellulose phase. It could be `Iβ`, `Iα`, `II` or `III`.

### Examples

```jldoctest

julia > atomsvecString(atoms, 4, 5, 7)
julia > atomsvecString("Iβ", 5, 7)
julia > atomsvecString("II")

```

"""

function atomsvecString(atomnames::Vector{String}, n::Int64, nblock::Int64)
    labels = Vector{String}[]
    atoms = repeat(replace.(atomnames, r"\d" => ""), outer=n);
    for blocks in collect(1:1:nblock)
        push!(labels, atoms)
    end
    return labels
end

function atomsvecString(phase::String, xsize::Int64, ysize::Int64)
    atomnames = map(_atomsvec_pick_first, get_crystallographic_info(phase)[1])
    if phase == "Ib" || phase == "Iβ"; n = 4; nblock=xsize*ysize; end
    if phase == "Ia" || phase == "Iα"; n = 2; nblock=2*xsize*ysize; end
    if phase == "II"; n = 4; nblock=xsize*ysize; end
    if phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"; n = 2; nblock=xsize*ysize; end
    return atomsvecString(atomnames, n, nblock), get_crystallographic_info(phase)[1]
end

function atomsvecString(phase::String)
    atomnames = map(_atomsvec_pick_first, get_crystallographic_info(phase)[1])
    if phase == "Ib" || phase == "Iβ"; n = 4; end
    if phase == "Ia" || phase == "Iα"; n = 2; end
    if phase == "II"; n = 4; end
    if phase == "III" || phase == "III_I" || phase == "III_i" || phase == "IIIi"; n = 2; end
    labels = repeat(replace.(atomnames, r"\d" => ""), outer=n);
    return labels, get_crystallographic_info(phase)[1]
end



"""
    _atomsvec_remove_numbers(string::String): removing the numbers inside the string::String.

"""
function _atomsvec_remove_numbers(string::String)
    return replace(string, r"\d" => "")
end

"""
    _atomsvec_pick_first(string::String): picking the first character in the string::String.

"""
function _atomsvec_pick_first(string::String)
    return "$(isempty(string) ? "" : string[1])"
end