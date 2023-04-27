import benchmarkedBurgers as burgers

burgers.main()

totalTime = 0
nTests = 127
for i in range(127):
    totalTime += burgers.main()
avg = totalTime / nTests
avgms = avg * 1000
print("Average time: " + str(avg))
print("In Milliseconds: " + str(avgms))
