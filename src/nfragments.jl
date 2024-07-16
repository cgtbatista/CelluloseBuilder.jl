function picking_fragments(vmdoutput)
    nfrag = nothing
    for id in eachindex(vmdoutput)
        str = string(vmdoutput[id])
        if occursin("Fragments:", str)
            nfrag = parse(Int64, split(str)[3])
        elseif isnothing(nfrag) && id == length(vmdoutput)
            error("The VMD output file does not contain the number of fragments.")
        end
    end
    return nfrag
end