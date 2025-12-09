const psfgen = joinpath(@__DIR__, "3rd-psfgen/psfgen")
const third_party_not_found = """
VMD/PSFGEN executable not found. Please, check if VMD/PSFGEN is installed or in your PATH.

Visit https://www.ks.uiuc.edu/Research/vmd/ for more information about VMD installation.
Visit https://www.ks.uiuc.edu/Research/psfgen/ for more information about PSFGEN installation.

Otherwise, as example, you can set the VMD executable path by running:
    ENV["VMD_EXECUTABLE"] = "/path/to/vmd"
"""

"""
    hasExecutable(exec::String)

Check if the executable is in the PATH and if it is executable.
"""
function hasExecutable(exec::String)
    path = Sys.which(exec)
    return isnothing(path) ? false : Sys.isexecutable(path)
end

"""
    execVMD(exec::String, input::String)

Execute VMD with the given input file and arguments. The output is returned as a string.
But for this julia package, the PSFGEN plugin is required and enough. So you can use the `execPSFGEN()` executable.
"""
function execVMD(exec::String, input::String)
    hasExecutable(exec) || error(third_party_not_found)
    isfile(input) || error("File not found: $input. Be sure to provide a readable file too.")
    try
        cmd = `$(exec) -dispdev text -e $input`
        return read(cmd, String)
    catch e
        throw(ErrorException("Error while executing VMD: $e")) # maybe sprint(showerror, e) can be prettier
    end
end

"""
    execPSFGEN(exec::String, input::String, args...)

Execute PSFGEN with the given input file and arguments. The output is returned as a string.
"""
function execPSFGEN(exec::String, input::String)
    hasExecutable(exec) || error(third_party_not_found)
    isfile(input) || error("File not found: $input. Be sure to provide a readable file too.")
    try
        cmd = `$(exec) $input`
        return read(cmd, String)
    catch e
        throw(ErrorException("Error while executing PSFGEN: $e")) # maybe sprint(showerror, e) can be prettier
    end
end

"""
    execPSFGEN(input::String)

Execute PSFGEN with the given input file using the default PSFGEN executable of this julia module.
"""
function execPSFGEN(input::String)
    return execPSFGEN(
        chmod(psfgen, 0o755),
        input
    )
end