include("DiscSimulations.jl")
using .DiscSimulations

initial_condition_sine(x, t, equation) = sinpi(x[1])
#at the moment, the burgers trixi sim is expecting a function that only takes x - this could be modified for more flexibility
initial_condition_burgers(x) = sinpi(x[1])

t_span = (0.0, 2.0)

source_zero(u, x, t, equations) = SVector(0, 0, 0)
params = DiscSimulations.Parameters(512, 0.0, 2π, initial_condition_burgers, source_zero, t_span)

solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
s = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

#NOTE: SHOW FUNCTION STILL NOT QUITE WORKING
#DiscSimulations.Base.show(IO, solution)

DiscSimulations.plotgif(s, t_span[1], t_span[2]; ylims = (-1, 1.1), legend = :topright)