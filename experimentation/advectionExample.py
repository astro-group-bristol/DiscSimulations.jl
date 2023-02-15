import math
import matplotlib.pyplot as pyplot

x0 = 0
t0 = 0
dx = 0.1
n = 4 * math.pi / dx
dt = 1
u = 0.5
data = {}

indices = []
x = x0

while x < x0 + n * dx:
    indices.append(x)
    x += dx

t = t0
data[t0] = [math.sin(x) for x in indices]

def step():
    global t
    prev = data[t]
    t = t + dt
    vals = []
    for i in range(len(indices)):
        x0 = prev[i]
        if i < len(indices) - 1:
            x1 = prev[i+1]
        else:
            # this is not the correct way to deal with a boundary, but
            # it works to demonstrate the principle
            x1 = prev[i]
        vals.append(x0 + (-1)*u*dt/dx*(x1-x0))
    data[t] = vals
 
step()
step()
pyplot.plot(data[t0])
pyplot.plot(data[t0 + dt])
pyplot.plot(data[t0 + 2 * dt])
pyplot.show()
