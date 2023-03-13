#modules and packages
include("./simpleGravity.jl")
using .simpleGravity
using Trixi
using OrdinaryDiffEq
using Plots
using Printf


function build_initial_conditions(user_supplied_ρ::Function)
    # define named lambda
    function _init_conditions(x, t, equations)
        # call user supplied function
        ρ = user_supplied_ρ(x, t, equations)
        #using the acceleration calculation in the og paper's setup for now
        q1 = -2*π*sin(π * x[1])
        return SVector(ρ,q1)
    end
    # closures are captured, so return the function to use in trixi
    return _init_conditions
end

#set up solver for equations
function setSolver(initCond, equations, mesh)
    boundary_conditions = (x_neg=BoundaryConditionDirichlet(((x, t, equations) -> SVector(sin(x[1]), -2*π*sin(π * x[1])))),
                       x_pos=BoundaryConditionDirichlet(((x, t, equations) -> SVector(sin(x[1]), -2*π*sin(π * x[1])))))
    solver = DGSEM(3, flux_central)
    semi = SemidiscretizationHyperbolic(mesh, equations, initCond, solver, boundary_conditions=boundary_conditions)
    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)
    return ode
end

#run simulations and collect results
function runAnimation(ode)
    callbacks = CallbackSet(SummaryCallback())
    sol = @time solve(ode, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7, callback=callbacks)
    tspan = collect(range(0.0, 1.0, 60))
    frames = Plots.@animate for t in tspan
        p = plot(sol(t))
        plot!(p, legend=:outertopright, title = Printf.@sprintf("t = %1.2f", t), xlims=(-80,80), xlabel = "x", ylabel="f")
    end
    gif(frames, "temp.gif", fps=10) |> display
end

function main(user_supplied_ρ::Function)
    #setting up grav eq
    equations = simpleGravity.gravEq1D()

    #arbitrary initial conditions
    x0 = 0
    t0 = 0
    mesh = TreeMesh(-1.0, 1.0, # min/max coordinates
                initial_refinement_level=4,
                periodicity=false,
                n_cells_max=10^4)

    #get initial conditions function and run
    _init_conditions(x, t, equations) = build_initial_conditions(user_supplied_ρ::Function)

    #solver
    ode = setSolver(_init_conditions(x0, t0, equations), equations, mesh)
    
    #run
    runAnimation(ode)
end

#from 3.1.1 of the paper. doesn't matter but thought i'd use it
#initial_condition_sine(x, t, equation::simpleGravity.gravEq1D) = 2 + (1/10) * sin(π*(x[1] - t))
initial_condition_sine(x, t, equation::simpleGravity.gravEq1D) = sin(x[1])

#other initial condition options:
#initial_condition = (x, t, equations) -> SVector(0.0)
#initial_condition_cos(x, t, equation::simpleGravity.gravEq1D) = SVector(2 + 2*cos(π*x[1]))
#initial_condition_x(x, t, equation::simpleGravity.gravEq1D) = x[1]

#main(initial_condition)

main(initial_condition_sine)