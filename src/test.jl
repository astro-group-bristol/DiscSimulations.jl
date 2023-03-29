include("DiscSimulations.jl")
using .DiscSimulationsjl

initial_condition_sine(x, t, equation) = sinpi(x[1])
source_zero(u, x, t, equations) = SVector(0, 0, 0)
params = DiscSimulationsjl.Parameters(512, 0, 2π, initial_condition_sine, source_zero)

solution = DiscSimulationsjl.main(params, DiscSimulationsjl.BurgersSimulation())
DiscSimulationsjl.plotgif(solution, 0.0, 2.0)