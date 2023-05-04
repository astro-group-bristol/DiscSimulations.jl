module DiscSimulations

using LinearAlgebra
using Plots
using Printf
using TruncatedStacktraces
using DifferentialEquations
TruncatedStacktraces.VERBOSE[] = true

#struct Parameters{T,DensityFuncType,SourceFuncType}
struct Parameters{T}
    xmin::T
    xmax::T
    ρ
    source
    t_span
    refinement_level::Int
    polydegree
    flux_type
end

abstract type SimulationType end
main(p::Parameters, ::SimulationType) = error("Unknown simulation type: $(typeof(t))")

struct BurgersSimulation <: SimulationType end
struct EulerWithSourceSimulation <: SimulationType end
struct SimpleGravitySimulation <: SimulationType end
struct SimpleBurgersSimulation <: SimulationType end
struct SpectralBurgersSimulation <: SimulationType end

main(p::Parameters, ::BurgersSimulation) = BurgerTrixi.solve_disc(p.xmin, p.xmax, p.refinement_level, p.ρ, p.t_span, p.polydegree, p.flux_type)
main(p::Parameters, ::EulerWithSourceSimulation) = EulerSource.solve_euler(p.source)
main(p::Parameters, ::SimpleGravitySimulation) = SimpleGravitySimulation.solve_gravity_sim(p.ρ)

function main(p::Parameters, ::SimpleBurgersSimulation)
    simple_x, _, simple_problem = BurgerSimple.setup(2^(p.refinement_level + 2), p.xmax, p.ρ, p.t_span)
    simple_sol = @time solve(simple_problem, Rodas5(autodiff = false); reltol = 1e-7, abstol = 1e-7)
    return simple_sol
end
function main(p::Parameters, ::SpectralBurgersSimulation)
    spectral_x, _, T, Ti, spectral_problem = BurgerSpectral.setup(2^(p.refinement_level + 2), p.xmax, p.ρ, p.t_span)
    spectral_sol = @time solve(spectral_problem, Rodas5(autodiff = false); reltol = 1e-7, abstol = 1e-7)
    return spectral_sol
end

include("solution.jl")

include("burgers/utils.jl")
include("burgers/trixi.jl")
include("burgers/spectral.jl")
include("burgers/simple.jl")
include("simpleGravitySim/simpleGravity.jl")
include("simpleGravitySim/simulation.jl")
#include("younsi-2012.jl")
include("Euler_with_source.jl")


export Parameters, SimulationType, BurgersSimulation, EulerWithSourceSimulation, SimpleGravitySimulation

end # module DiscSimulations
