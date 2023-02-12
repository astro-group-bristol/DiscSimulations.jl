module BurgerSimple

using DifferentialEquations

import DiscSimulations: STANDARD_BURGER_INIT, STANDARD_BURGER_TSPAN, central_difference

function burger(du, u, p, t)
    # unpack additional paramters
    D = p
    # dirichlet boundary condition
    u[1] = u[end] = 0.0
    step = - u .* (D * u)
    du .= step
end

function setup(N, x_max, init = STANDARD_BURGER_INIT, t_span = STANDARD_BURGER_TSPAN)
    Δx = x_max / N

    x = collect(range(0.0, x_max, N))
    u0 = init(x) 

    # get central difference operator as matrix
    D = central_difference(Δx, N)

    problem = ODEProblem{true}(
        burger,
        u0,
        t_span,
        D
    )

    x, u0, problem
end

end # module BurgerSimple

export BurgerSimple