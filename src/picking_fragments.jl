"""

    picking_fragments(vmdoutput)

Picking the number of fragments (e.g. cellulose chains) inside the structure loaded on VMD. Usualy,
the VMD output has the type `Vector{SubString{String}}` such as:

```julia
1817-element Vector{SubString{String}}:
 "Info) VMD for LINUXAMD64, version 1.9.4a57 (April 27, 2022)"
 "Info) http://www.ks.uiuc.edu/Research/vmd/                         "
 "Info) Email questions and bug reports to vmd@ks.uiuc.edu           "
 "Info) Please include this reference in published work using VMD:   "
 "Info)    Humphrey, W., Dalke, A. and Schulten, K., `VMD - Visual   "
 "Info)    Molecular Dynamics', J. Molec. Graphics 1996, 14.1, 33-38."
 "Info) -------------------------------------------------------------"
 "Info) Multithreading available, 24 CPUs."
 "Info)   CPU features: SSE2 SSE4.1 AVX AVX2 FMA F16 HT "
 "Info) Free system memory: 19GB (61%)"
 "Info) Creating CUDA device pool and initializing hardware..."
 "Info) Detected 1 available CUDA accelerator::"
 "Info) [0] NVIDIA GeForce RTX 40" ⋯ 18 bytes ⋯ "5 GHz, 7.7GB RAM SP64 KT AE2 ZC"
 "Info) Detected 1 available TachyonL/OptiX ray tracing accelerator"
 "Info)   Compiling  OptiX shaders on 1 target GPU..."
 "Info) Dynamically loaded 3 plugins in directory:"
 "Info) /usr/local/lib/vmd/plugins/LINUXAMD64/molfile"
 "2.0"
 "psfgen) reading topology file /" ⋯ 47 bytes ⋯ "Z/src/toppar/top_all36_carb.rtf"
 ""
 "psfgen)  \$Id: top_allxx_sugar.inp,v 1.106 2014/08/19 19:07:43 alex Exp \$"
 "psfgen) >>>>>>>>>>>> All-hydrogen topologies used in the <<<<<<<<<<<<<<<<"
 "psfgen) >>>>> development of the CHARMM carbohydrate  force field<<<<<<<<"
 ""
 ⋮
 "Info) VMD for LINUXAMD64, version 1.9.4a57 (April 27, 2022)"
 "Info) Exiting normally."
 ""
```
## Arguments

- `vmdoutput::Vector{SubString{String}}`: output of the structure loaded on VMD.

### Examples

```jldoctest

julia > picking_fragments(vmdoutput)
julia > 59

```

"""
function picking_fragments(vmdoutput)

    nfragments = nothing
    
    for id in eachindex(vmdoutput)
        str = string(vmdoutput[id])
        if occursin("Fragments:", str)
            nfragments = parse(Int64, split(str)[3])
        elseif isnothing(nfragments) && id == length(vmdoutput)
            error("The VMD output file does not contain the number of fragments.")
        end
    end
    
    return nfragments
    
end