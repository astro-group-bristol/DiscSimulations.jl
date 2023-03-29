module simpleGravitySimulation

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
        #for ρ_v1 and ρ_e, using arbitrary values for now
        ρ_v1 = ρ * 0.1
        ρ_e = ρ * 10
        a1 = 0.0
        SVector(ρ, ρ_v1, ρ_e, a1)
    end
    # closures are captured, so return the function to use in trixi
    return _init_conditions
end

#set up solver for equations
function setSolver(initCond, equations, mesh)
    boundary_conditions = (x_neg=BoundaryConditionDirichlet(Trixi.boundary_condition_do_nothing),
                       x_pos=BoundaryConditionDirichlet(Trixi.boundary_condition_do_nothing))
    solver = DGSEM(3, flux_hll)
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

function solve_gravity_sim(user_supplied_ρ::Function)
    γ = 2.0
    #setting up grav eq
    equations = simpleGravity.GravityEquations1D(γ)
    #get initial conditions function and run
    _init_conditions = build_initial_conditions(user_supplied_ρ::Function)
    mesh = TreeMesh(XMIN, XMAX, # min/max coordinates
                initial_refinement_level=4,
                periodicity=false,
                n_cells_max=N)
    ode, semi = setSolver(_init_conditions, equations, mesh)
    sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7, callback=callbacks)
    DiscSolution(sol, semi, OneDimension())
end

end #module

export simpleGravitySimulation