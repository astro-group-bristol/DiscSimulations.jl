using Plots

x0 = float(0)
t0 = 0
dx = 0.1
n = 4 * π / dx
dt = 1
u = 0.5
data = Dict()

indices = []
x = x0

while x < x0 + n * dx
    push!(indices, x)
    x += dx
end

t = t0
temp = []
for x in indices
    push!(temp, sin(x))
end
data[t0] = temp

function step()
    global t
    prev = data[t]
    t = t + dt
    vals = []
    for i in 1:length(indices)
        x0 = prev[i]
        if i < length(indices) - 1
            x1 = prev[i + 1]
        else
            x1 = prev[i]
        end
        push!(vals, (x0 + (-1)*u*dt/dx*(x1 - x0)))
    end
    data[t] = vals
end

step()
step()
plot(data[t0])
plot!(data[t0 + dt])
plot!(data[t0 + 2*dt])