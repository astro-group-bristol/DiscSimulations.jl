module BurgerTrixi

using Trixi
using DifferentialEquations
using Plots
using Printf

include("utils.jl")
using .BurgerUtils

include("../solution.jl")

#import DiscSimulations:
#    BurgerUtils.STANDARD_BURGER_INIT, BurgerUtils.STANDARD_BURGER_TSPAN, DiscSolution, OneDimension

function setup(N, x_max, init = BurgerUtils.STANDARD_BURGER_INIT, t_span = BurgerUtils.STANDARD_BURGER_TSPAN)
    x_min = 0.0

    equations = InviscidBurgersEquation1D()
    solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

    mesh = TreeMesh(x_min, x_max, initial_refinement_level = 3, n_cells_max = N)
    semi = SemidiscretizationHyperbolic(mesh, equations, (x, t, eqs) -> init(x), solver)

    problem = semidiscretize(semi, t_span)

    semi, problem
end

function solve_disc(N, x_max, init = BurgerUtils.STANDARD_BURGER_INIT, t_span = BurgerUtils.STANDARD_BURGER_TSPAN)
    semi, problem = setup(N, x_max, init, t_span)
    sol = solve(problem, RDPK3SpFSAL49(), abstol = 1.0e-7, reltol = 1.0e-7)
    return DiscSolution(sol, semi, OneDimension())
end

end # module BurgerTrixi

export BurgerTrixi
