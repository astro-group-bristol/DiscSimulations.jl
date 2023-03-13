module BurgerSpectral

using ApproxFun
using LinearAlgebra
using DifferentialEquations

import DiscSimulations: STANDARD_BURGER_INIT, STANDARD_BURGER_TSPAN

function burger_fourier(du_f, u_f, p, t)
    # unpack parameters
    D, T, Ti, tmp1, tmp2 = p
    # calculate derivative in fourier space
    mul!(tmp1, D, u_f)
    # transform derivative to regular space
    mul!(tmp2, Ti, tmp1)
    # transform u to regular space
    mul!(tmp1, Ti, u_f)
    # D * u
    @. tmp1 = tmp2 * tmp1
    # transform back to fourier space
    mul!(tmp2, T, tmp1)
    du_f .= -tmp2
end

function setup(N, x_max, init = STANDARD_BURGER_INIT, t_span = STANDARD_BURGER_TSPAN)
    space = Fourier()
    x = points(space, N)

    D = (Derivative(space) → space)[1:N, 1:N]

    # fourier transform
    T = ApproxFun.plan_transform(space, N)
    # inverse fourier transform
    Ti = ApproxFun.plan_itransform(space, N)

    u0 = T * init(x)
    # pre-allocate temporary arrays
    tmp1 = similar(u0)
    tmp2 = similar(u0)

    problem = ODEProblem{true}(burger_fourier, u0, t_span, (D, T, Ti, tmp1, tmp2))

    x, u0, T, Ti, problem
end

end # module BurgerSpectral

export BurgerSpectral