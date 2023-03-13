module BurgerTrixi

using Trixi

import DiscSimulations: STANDARD_BURGER_INIT, STANDARD_BURGER_TSPAN

function setup(N, x_max, init = STANDARD_BURGER_INIT, t_span = STANDARD_BURGER_TSPAN)
    x_min = 0.0

    equations = InviscidBurgersEquation1D()
    solver = DGSEM(polydeg=3, surface_flux = flux_lax_friedrichs)

    mesh = TreeMesh(x_min, x_max, initial_refinement_level=3, n_cells_max = N)
    semi = SemidiscretizationHyperbolic(mesh, equations, (x, t, eqs) -> init(x), solver)

    problem = semidiscretize(semi, t_span)
    
    semi, problem
end

end # module BurgerTrixi

export BurgerTrixi