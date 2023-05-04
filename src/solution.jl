import Trixi: PlotData1D, PlotData2D

abstract type AbstractDimension end
struct TwoDimension <: AbstractDimension end
struct OneDimension <: AbstractDimension end

struct DiscSolution{Dimension,SolutionType,DiscretiziationType}
    sol::SolutionType
    semi::DiscretiziationType
end

function DiscSolution(sol, semi, dimension_type::AbstractDimension)
    DiscSolution{typeof(dimension_type),typeof(sol),typeof(semi)}(sol, semi)
end

#1D version of plot
function dataplot(solution::DiscSolution{<:OneDimension}, solver, t)
    return PlotData1D(solver(t), solution.semi)
end

#2D version of plot
function dataplot(solution::DiscSolution{<:TwoDimension}, solver, t)
    return PlotData2D(solver(t), solution.semi)
end

dataplot(solution::DiscSolution{<:AbstractDimension}, t) = error("Unknown simulation type: $(typeof(t))")

function plotgif(solution::DiscSolution, tmin, tmax; kwargs...)
    ts = range(tmin, tmax, 150)
    solver = @time solution.sol
    frames = Plots.@animate for t in ts
        pd = dataplot(solution, solver, t)
        plot(pd ; kwargs...)
    end
    gif(frames, "temp.gif", fps = 10) |> display
end

function Base.show(io::IO, sol::DiscSolution, x)
    println(io, string(typeof(sol).parameters[x]))
end