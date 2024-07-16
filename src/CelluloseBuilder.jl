module CelluloseBuilder

export cellulosebuilder

# Write your package code here.
const DEFAULT_CARB_TOPOLOGY_FILE = "$(@__DIR__)/toppar/top_all36_carb.rtf"

# main cellulose builder function
include("./cellulose.jl")

include("./atomnames.jl")
include("./cleaning.jl")
include("./crystallographic.jl")
include("./formatter.jl")
include("./nfragments.jl")
include("./pbc.jl")
include("./systemgenerator.jl")

end