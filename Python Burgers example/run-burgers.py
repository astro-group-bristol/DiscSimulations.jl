import benchmarkedBurgers as burgers
import math

totalTime = 0
nTests = 3
results = []
for i in range(nTests):
    time = burgers.main()
    results.append(time)
    totalTime += time
avg = totalTime / nTests
totaldif = 0
for i in results:
    totaldif += (i - avg)**2
sd = math.sqrt(totaldif / (totalTime - 1))

avgms = avg * 1000
print("Average time: " + str(avg))
print("In Milliseconds: " + str(avgms))
print("SD: " + str(sd))
