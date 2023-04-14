module DiscSimulations

using LinearAlgebra
using Plots
using Printf
using TruncatedStacktraces
TruncatedStacktraces.VERBOSE[] = true

#struct Parameters{T,DensityFuncType,SourceFuncType}
struct Parameters{T}
    N::Int
    xmin::T
    xmax::T
    ρ
    source
end

abstract type SimulationType end
main(p::Parameters, ::SimulationType) = error("Unknown simulation type: $(typeof(t))")

struct BurgersSimulation <: SimulationType end
struct EulerWithSourceSimulation <: SimulationType end
struct SimpleGravitySimulation <: SimulationType end

main(p::Parameters, ::BurgersSimulation) = BurgerTrixi.solve_disc(p.N, p.xmax)
main(p::Parameters, ::EulerWithSourceSimulation) = EulerSource.solve_euler(p.source)
main(p::Parameters, ::SimpleGravitySimulation) = SimpleGravitySimulation.solve_gravity_sim(p.ρ)

include("solution.jl")

include("burgers/utils.jl")
include("burgers/trixi.jl")
#include("burgers/spectral.jl")
#include("burgers/simple.jl")
include("simpleGravitySim/simpleGravity.jl")
include("simpleGravitySim/simulation.jl")
include("younsi-2012.jl")
include("Euler_with_source.jl")


export Parameters, SimulationType, BurgersSimulation, EulerWithSourceSimulation, SimpleGravitySimulation

end # module DiscSimulations
