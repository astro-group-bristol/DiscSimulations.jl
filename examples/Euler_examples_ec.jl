using OrdinaryDiffEq
using Trixi
using Plots


###############################################################################
# setting up source term
#function of form S(Q) where Q = (ρ, ρux, ρuy, ρE)
@inline function source_test(u, x, t, equations)
    #in 1D for now so we are leaving out y
    ϕ = tan(x[1])^(-1) #arctan(y/x)
    ρ = u[1] #density
    G = 6.67408 * 10^(-11) #gravitational constant
    mₛ = 4 * 1.989 * 10^(30) #mass of central object
    ux = u[2] #velocity of x
    uy = u[3] #velocity of y
    r = sqrt(x[1]^2) #√(x² + y²)

    du1 = -cos(ϕ)*ρ*(G*mₛ/(r^2))
    du2 = -sin(ϕ)*ρ*(G*mₛ/(r^2))
    du3 = -(ux*cos(ϕ) + uy*sin(ϕ))*ρ*(G*mₛ/(r^2))
    return SVector(du1, du2, du3)
end

source_zero(u, x, t, equations) = SVector(0, 0, 0)

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_weak_blast_wave

volume_flux = flux_ranocha
solver = DGSEM(polydeg=3, surface_flux=flux_ranocha,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux), source_terms=source_zero)

coordinates_min = (-2.0,)
coordinates_max = ( 2.0,)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=5,
                n_cells_max=10_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.8)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

plot(sol)