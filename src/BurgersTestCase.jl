include("DiscSimulations.jl")
using .DiscSimulations
using BenchmarkTools

xmin = 0.0
xmax = 1.0
N = 1024

tmax = (xmax - xmin)/1.0
t_span = (0.0, 0.2*tmax)

function initial_condition_burgers(x)
    u = 1.0
    if(x[1] >= 0.333 && x[1] <= 0.666)
        u = u + 0.5*sin(2.0*π*(x[1]-0.333)/0.333)
    end
    return u
end

#params expects this so we're defining a dummy one
source_zero(u, x, t, equations) = SVector(0, 0, 0)

params = DiscSimulations.Parameters(N, xmin, xmax, initial_condition_burgers, source_zero, t_span, 8)

@benchmark solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
s = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

open("julia-output.csv", "w") do file
    write(file, string(last(s.sol)))
end

DiscSimulations.plotgif(s, t_span[1], t_span[2]; xlims = (0, 1), ylims = (0, 2), legend = :topright)