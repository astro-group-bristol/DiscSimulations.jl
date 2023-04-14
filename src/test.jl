include("DiscSimulations.jl")
using .DiscSimulations

initial_condition_sine(x, t, equation) = sinpi(x[1])
source_zero(u, x, t, equations) = SVector(0, 0, 0)
params = DiscSimulations.Parameters(512, 0.0, 2π, initial_condition_sine, source_zero)

solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
s = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
kwargs = ["ylims = (-1, 1.1)", "legend = :topright"]
DiscSimulations.plotgif(s, 0.0, 2.0)#, kwargs)