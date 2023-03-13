module DiscSimulations

using LinearAlgebra

include("burgers/utils.jl")
include("burgers/spectral.jl")
include("burgers/simple.jl")
include("burgers/trixi.jl")
include("younsi-2012.jl")

end # module DiscSimulations
