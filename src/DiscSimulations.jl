module DiscSimulations

using LinearAlgebra
using Plots
using Printf

include("solution.jl")

include("burgers/utils.jl")
include("burgers/spectral.jl")
include("burgers/simple.jl")
include("burgers/trixi.jl")
include("younsi-2012.jl")


end # module DiscSimulations
