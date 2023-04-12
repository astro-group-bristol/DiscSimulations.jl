include("src/DiscSimulations.jl")
using .DiscSimulations

initial_condition_sine(x, t, equation) = sinpi(x[1])
source_zero(u, x, t, equations) = SVector(0, 0, 0)
params = DiscSimulations.Parameters(512, 0, 2π, initial_condition_sine, source_zero)

solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
DiscSimulations.plotgif(solution, 0.0, 2.0)