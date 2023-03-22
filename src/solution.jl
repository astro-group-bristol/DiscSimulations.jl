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

function plotgif(sol::DiscSolution{OneDimension}, tmin, tmax)
    ts = range(tmin, tmax, 150)
    frames = Plots.@animate for t in ts
        # 1d
        pd = PlotData1D(sol.sol(t), sol.semi)
        plot(
            pd,
            label = Printf.@sprintf("t = %1.2f", t),
            ylims = (-0.1, 2.1),
            legend = :topright,
        )
    end
    gif(frames, "temp.gif", fps = 10) |> display
end

function plotgif(sol::DiscSolution{TwoDimension}, tmin, tmax)
    ts = range(tmin, tmax, 150)
    frames = Plots.@animate for t in ts
        # 2d
        pd = PlotData2D(sol.sol(t), sol.semi)
        plot(
            pd,
            label = Printf.@sprintf("t = %1.2f", t),
            ylims = (-0.1, 2.1),
            legend = :topright,
        )
    end
    gif(frames, "temp.gif", fps = 10) |> display
end
