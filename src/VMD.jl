function execVMD(vmd::String, input::String)
    
    if !hasVMD(vmd)
        error("VMD executable not found. Please, check if VMD is installed and in your PATH.")
    end

    if !isfile(input)
        error("File not found: $input. Be sure to provide a valid file.")
    end

    return Base.read(`$vmd -dispdev text -e $input`, String)
end

function hasVMD(vmd::String)
    return Sys.which(vmd) !== nothing
end