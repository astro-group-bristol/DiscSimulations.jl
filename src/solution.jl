import Trixi: PlotData1D, PlotData2D

abstract type AbstractDimension end
struct TwoDimension <: AbstractDimension end
struct OneDimension <: AbstractDimension end

struct DiscSolution{Dimension,SolutionType,DiscretiziationType}
    sol::SolutionType
    semi::DiscretiziationType
    #to add: keyword args
end

function DiscSolution(sol, semi, dimension_type::AbstractDimension)
    DiscSolution{typeof(dimension_type),typeof(sol),typeof(semi)}(sol, semi)
end

#1D version of plot
function dataplot(solution::DiscSolution{OneDimension})
    return PlotData1D(solution.sol(t), solution.semi)
end

#2D version of plot
function dataplot(solution::DiscSolution{TwoDimension})
    return PlotData2D(solution.sol(t), solution.semi)
end

function plotgif(solution, tmin, tmax)
    ts = range(tmin, tmax, 150)
    frames = Plots.@animate for t in ts
        pd = dataplot(solution)
        plot(
            pd,
            label = Printf.@sprintf("t = %1.2f", t),
            ylims = (-0.1, 2.1),
            legend = :topright,
        )
    end
    gif(frames, "temp.gif", fps = 10) |> display
end
