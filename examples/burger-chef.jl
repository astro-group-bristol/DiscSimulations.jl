# welcome to burger king can i take your order please
# solving Burger's equations N ways
using DiscSimulations
using Plots
using DifferentialEquations
using Printf

N = 512
xmax = 2π

simple_x, _, simple_problem = BurgerSimple.setup(N, xmax)
simple_sol = @time solve(simple_problem, Rodas5(autodiff=false); reltol = 1e-7, abstol = 1e-7)

spectral_x, _, T, Ti, spectral_problem = BurgerSpectral.setup(N, xmax)
spectral_sol = @time solve(spectral_problem, Rodas5(autodiff=false); reltol = 1e-7, abstol = 1e-7)

# trixi_x, trixi_problem = BurgerTrixi.setup(N, xmax)
# trixi_sol = @time solve(trixi_problem, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7)

# use the new disc solution type
discsol = BurgerTrixi.solve_disc(N, xmax)
# plot that individually
DiscSimulations.plotgif(discsol, 0.0, 1.0)

tspan = collect(range(0.0, 1.0, 60))
frames = Plots.@animate for t in tspan
    p = plot()

    # plot!(trixi_sol(t), trixi_x, label="Trixi")
    plot!(simple_x, simple_sol(t), label="Simple")
    plot!(spectral_x, Ti * spectral_sol(t), label="Spectral")

    plot!(p, legend=:outertopright, title = Printf.@sprintf("t = %1.2f", t), ylims = (-0.1, 2.1), xlabel = "x", ylabel="f")
end

gif(frames, "temp.gif", fps=10) |> display