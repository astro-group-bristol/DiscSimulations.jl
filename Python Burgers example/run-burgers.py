import benchmarkedBurgers as burgers

totalTime = 0
nTests = 1
for i in range(nTests):
    totalTime += burgers.main()
avg = totalTime / nTests
avgms = avg * 1000
print("Average time: " + str(avg))
print("In Milliseconds: " + str(avgms))


