#modules and packages
include("./simpleGravity.jl")
using .simpleGravity
using Trixi
using OrdinaryDiffEq
using Plots
using Printf

T_SPAN = 5.0
XMIN = 0
XMAX = 2π
N = 512

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
    boundary_conditions = (x_neg=BoundaryConditionDirichlet(Trixi.boundary_condition_do_nothing),
                       x_pos=BoundaryConditionDirichlet(Trixi.boundary_condition_do_nothing))
    solver = DGSEM(3, flux_central)
    semi = SemidiscretizationHyperbolic(mesh, equations, initCond, solver, boundary_conditions=boundary_conditions)
    tspan = (0.0, T_SPAN)
    ode = semidiscretize(semi, tspan)
    return ode, semi
end

function runAnimation(ode, semi)
    callbacks = CallbackSet(SummaryCallback())
    sol = @time solve(ode, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7, callback=callbacks)
    begin
        ts = collect(range(0.0, T_SPAN, 150))
        frames = Plots.@animate for t in ts
            pd = PlotData1D(sol(t), semi)
            # plot!(p, legend=:outertopright, title = Printf.@sprintf("t = %1.2f", t), xlabel = "x", ylabel="f", ylims = (-1, 1))
            plot(pd, label = Printf.@sprintf("t = %1.2f", t), ylims = (-1, 1.1), legend = :topright)
        end
        gif(frames, "temp.gif", fps=10) |> display
    end
end

function main(user_supplied_ρ::Function)
    #setting up grav eq
    equations = simpleGravity.gravEq1D()

    mesh = TreeMesh(XMIN, XMAX, # min/max coordinates
                initial_refinement_level=4,
                periodicity=false,
                n_cells_max=N)

    #get initial conditions function and run
    _init_conditions = build_initial_conditions(user_supplied_ρ::Function)
    #solver
    ode, semi = setSolver(_init_conditions, equations, mesh)
    
    #run
    runAnimation(ode, semi)
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