#modules and packages
using .simpleGravity
using Trixi
using OrdinaryDiffEq
using Plots
using Printf

function get_accelerations(x, t, equations)
    #original eqs factored in y but that is 0 here
    q1 = -2*π*sin(π*x)*sin(0)
    q2 = 4*π*cos(π*x)*cos(0)
    return q1, q2
end

function build_initial_conditions(user_supplied_ρ::Function)
    # define named lambda
    function _init_conditions(x, t, equations)
        # call user supplied function
        ρ = user_supplied_ρ(x, t, equations)
        # acceleration fields
        acc_x, acc_y = get_accelerations(x, t, equations)
        return SVector(ρ, acc_x, acc_y)
    end
    # closures are captured, so return the function to use in trixi
    return _init_conditions
end

#set up solver for equations
function setSolver(initCond, mesh)
    solver = DGSEM(3, flux_central)
    semi = SemidiscretizationHyperbolic(mesh, equation, initial_condition_sine, solver)
    tspan = (0.0, 0.1)
    ode = semidiscretize(semi, tspan)
    return ode
end

#run simulations and collect results
function runAnimation(ode)
    sol = @time solve(ode, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7)
    tspan = collect(range(0.0, 1.0, 60))
    frames = Plots.@animate for t in tspan
        p = plot(sol(t))
        plot!(p, legend=:outertopright, title = Printf.@sprintf("t = %1.2f", t), xlabel = "x", ylabel="f")
    end
    gif(frames, "temp.gif", fps=10) |> display
end

function main(user_supplied_ρ::Function)
    #equations
    gravEq = simpleGravity.gravEq()
    x = ##
    t = ##
    equations = ##

    #get initial conditions function and run
    _init_conditions(x, t, equations) = build_initial_conditions(user_supplied_ρ::Function)
    initCond = _init_conditions(x, t, equations)

    #solver
    ode = setSolver(initCond, mesh)
    
    #run
    runAnimation(ode)
end

main(sin())