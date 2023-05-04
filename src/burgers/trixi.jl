module BurgerTrixi

using Trixi
using DifferentialEquations
using Plots
using Printf

import DiscSimulations:
    STANDARD_BURGER_INIT, STANDARD_BURGER_TSPAN, DiscSolution, OneDimension

function setup(x_min, x_max, init, t_span, refinement_level, polydegree, flux_type)
    N = 2^(refinement_level + 2)
    equations = InviscidBurgersEquation1D()
    solver = DGSEM(polydeg = polydegree, surface_flux = flux_type)

    mesh = TreeMesh(x_min, x_max, initial_refinement_level = refinement_level, n_cells_max = N)
    semi = SemidiscretizationHyperbolic(mesh, equations, (x, t, eqs) -> init(x), solver)

    ode = semidiscretize(semi, t_span)

    semi, ode
end

function solve_disc(x_min, x_max, refinement_level, init = STANDARD_BURGER_INIT, t_span = STANDARD_BURGER_TSPAN, polydegree = 3, flux_type = flux_lax_friedrichs)
    semi, problem = setup(x_min, x_max, init, t_span, refinement_level, polydegree, flux_type)
    sol = solve(problem, RDPK3SpFSAL49(), abstol = 1.0e-7, reltol = 1.0e-7)
    return DiscSolution(sol, semi, OneDimension())
end

end # module BurgerTrixi

export BurgerTrixi
