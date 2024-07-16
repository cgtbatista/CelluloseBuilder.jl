function cleaning_tmpfiles()
    files = readdir("/tmp")
    for file in files
        if occursin("jl_", file) || occursin("tmp_", file) || occursin(".tcl", file)
            rm("/tmp/$file")
        end
    end
    return nothing
end