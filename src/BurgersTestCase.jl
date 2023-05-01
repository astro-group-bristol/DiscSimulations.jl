include("DiscSimulations.jl")
using .DiscSimulations
using BenchmarkTools
using Trixi
using Plots

xmin = 0.0
xmax = 1.0

tmax = (xmax - xmin)/1.0
t_span = (0.0, 0.2*tmax)

function initial_condition_burgers(x)
    u = 1.0
    if(x[1] >= 0.333 && x[1] <= 0.666)
        u = u + 0.5*sin(2.0*π*(x[1]-0.333)/0.333)
    end
    return u
end

#params expects this so we're defining a dummy one
source_zero(u, x, t, equations) = SVector(0, 0, 0)

params = DiscSimulations.Parameters(xmin, xmax, initial_condition_burgers, source_zero, t_span, 10, 3, Trixi.flux_lax_friedrichs)

@benchmark solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
s = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

s1 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s2 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s3 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s4 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s5 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s6 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s7 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s8 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s9 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
s10 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())


fg = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
feo = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

poly2 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
poly4 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
poly5 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
poly6 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
poly7 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())
poly8 = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

function importAndRead(filexs, fileus)
    pyDXs = []
    pyDUs = []
 
    #import python data
    pyDataX = open(filexs, "r")
    line = 0  
    while ! eof(pyDataX)          
       s = parse(Float64, readline(pyDataX))         
       line += 1
       push!(pyDXs, s)
    end #while
    close(pyDataX)
 
    pyDataU = open(fileus, "r")
    line = 0  
    while ! eof(pyDataU)          
       s = parse(Float64, readline(pyDataU)) 
       line += 1
       push!(pyDUs, s)
    end #while
    close(pyDataU)
 
    return pyDXs, pyDUs
 end 

pyDXs1024, pyDUs1024 = importAndRead("python-xs-1024.txt", "python-us-1024.txt")
pyDXs256, pyDUs256 = importAndRead("python-xs-256.txt", "python-us-256.txt")

pd = PlotData1D(solution.sol)

plot(pd, title="Trixi and Python Burgers Comparison", ylabel="u")
plot!(pyDXs, pyDUs)
savefig("PythonVsTrixi-1024.png")

pd1 = PlotData1D(s1.sol)
pd2 = PlotData1D(s2.sol)
pd3 = PlotData1D(s3.sol)
pd4 = PlotData1D(s4.sol)
pd5 = PlotData1D(s5.sol)
pd6 = PlotData1D(s6.sol)
pd7 = PlotData1D(s7.sol)
pd8 = PlotData1D(s8.sol)
pd9 = PlotData1D(s9.sol)
pd10 = PlotData1D(s10.sol)

pdfg = PlotData1D(fg.sol)
pdfeo = PlotData1D(feo.sol)

pdpoly2 = PlotData1D(poly2.sol)
pdpoly4 = PlotData1D(poly4.sol)
pdpoly5 = PlotData1D(poly5.sol)
pdpoly6 = PlotData1D(poly6.sol)
pdpoly7 = PlotData1D(poly7.sol)
pdpoly8 = PlotData1D(poly8.sol)

plot(pdpoly8, title="Trixi Polynomial 8", ylabel="u")
plot!(pyDXs, pyDUs)
savefig("TrixiPoly8.png")
savefig("TrixiPoly8Comp.png")

plot(pd10, title="Trixi Initial Refinement 10", ylabel="u", label="Trixi")
plot!(pyDXs256, pyDUs256, label="Python")
plot!(legend=:topleft)
savefig("InitRefinement10Compare.png")

#error calculation
t_dy = 0
t_abs_dy = 0   
t_relerr = 0  
t_pererr = 0   
t_mean_err = 0    
t_MSE = 0        
t_RMSE = 0

y = last(s7.sol)
inc = length(y) / 256
inc = 256 / length(y)
for i in 1:256
    #y0 = pyDUs256[trunc(Int, i*inc)]
    #y1 = y[i]

    y0 = pyDUs256[i]
    y1 = y[trunc(Int, i*inc)]

    dy = y0-y1 # error 
    abs_dy = abs(y0-y1)   # absolute error 
    relerr = abs(y0-y1)./y0  # relative error 
    pererr = abs(y0-y1)./y0*100   # percentage error 
    mean_err = mean(abs(y0-y1))    # mean absolute error 
    MSE = mean((y0-y1).^2)        # Mean square error 
    #RMSE = sqrt(mean((y0-y1).^2)) # Root mean square error

    t_dy += dy
    t_abs_dy += abs_dy   
    t_relerr += relerr  
    t_pererr += pererr   
    t_mean_err += mean_err    
    t_MSE += MSE        
    #t_RMSE += RSME
end

function printErrs(t_dy, t_abs_dy, t_relerr, t_pererr, t_mean_err, t_MSE, N)
    print("dy: ", string(t_dy / N), "\n")
    print("abs_dy: ", string(t_abs_dy / N), "\n")   
    print("relerr: ", string(t_relerr / N), "\n")  
    print("pererr: ", string(t_pererr / N), "\n")   
    print("mean_err: ", string(t_mean_err / N), "\n")    
    print("MSE: ", string(t_MSE / N), "\n")        
end

printErrs(t_dy, t_abs_dy, t_relerr, t_pererr, t_mean_err, t_MSE, 256)

print(string(solver(tmax)))

open("julia-output.txt", "w") do file
    write(file, string(solution.sol))
end

DiscSimulations.plotgif(poly2, t_span[1], t_span[2]; xlims = (0, 1), ylims = (0, 2), legend = :topright)