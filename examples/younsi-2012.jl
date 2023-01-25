using DiscSimulations
using Plots

# black hole mass
M = 1.0

# black hole spin
a = 0.2

# keplerian radius
rₖ = 12.0

# angular velocity profile index of the torus
n = 0.21

begin
    p = plot(xlabel = "r", ylabel = "z", legend = :outertopright)
    for a in [0.0, 0.2, 0.4, 0.6, 0.8]
        r, z = @time isobar(M, a, n, rₖ)
        plot!(r, z, label = "a = $a")
        # plot but upside down
        plot!(r, -z, label = "a = $a")
    end
    p
end
