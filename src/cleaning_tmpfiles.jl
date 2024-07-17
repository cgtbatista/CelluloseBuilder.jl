"""

    cleaning_tmpfiles()

Move the exported files to the work directory and clean all the *.pdb, *.tcl, and *xyz on temporary
files inside the temp folder (it can be different on Linux, Windows and MacOS).

### Arguments

- `filename::String`: main name of the exported file.
- `destination_path`: the default destination folder is the working directory.

### Examples

```jldoctest

julia > cleaning_tmpfiles()

```

"""

function cleaning_tmpfiles(filename::String; destination_path=pwd())
    
    tmp_path = tempdir()

    mv(joinpath(tmp_path, "$filename.xyz"), joinpath(destination_path, "$filename.xyz"), force=true)
    mv(joinpath(tmp_path, "$filename.psf"), joinpath(destination_path, "$filename.psf"), force=true)
    mv(joinpath(tmp_path, "$filename.pdb"), joinpath(destination_path, "$filename.pdb"), force=true)

    files = readdir(tmp_path)
    for file in files
        if occursin(".pdb", file) || occursin(".xyz", file) || occursin(".tcl", file)
            rm(joinpath(tmp_path, file))
        end
    end

    return nothing

end