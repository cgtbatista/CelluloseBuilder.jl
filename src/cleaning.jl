function cleaning_tmpfiles()
    files = readdir("/tmp")
    for file in files
        if occursin(".pdb", file) || occursin(".xyz", file) || occursin(".tcl", file)
            rm("/tmp/$file")
        end
    end
    return nothing
end