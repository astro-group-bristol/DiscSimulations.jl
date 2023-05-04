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

params = DiscSimulations.Parameters(xmin, xmax, initial_condition_burgers, source_zero, t_span, 8, 3, Trixi.flux_lax_friedrichs)

@benchmark solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
solution = DiscSimulations.main(params, DiscSimulations.BurgersSimulation())
s = DiscSimulations.DiscSolution(solution.sol, solution.semi, DiscSimulations.OneDimension())

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

#error calculation
t_dy = 0
t_abs_dy = 0   
t_relerr = 0  
t_pererr = 0   
t_mean_err = 0    
t_MSE = 0        
t_RMSE = 0

y = last(poly4.sol)
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

DiscSimulations.plotgif(poly2, t_span[1], t_span[2]; xlims = (0, 1), ylims = (0, 2), legend = :topright)