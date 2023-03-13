module simpleGravity

using Trixi, OrdinaryDiffEq, LinearAlgebra, ForwardDiff;

struct gravEq1D <: Trixi.AbstractEquations{1 #= number of spatial dimensions =#,
                                        2 #= number of primary variables, (Φ, q1)=#}
end

@inline Trixi.flux(u, orientation, equation::gravEq1D) = SVector(0, -u[2])
#@inline Trixi.flux(u, orientation, equation::gravEq1D) = u
Trixi.varnames(::typeof(cons2cons), ::gravEq1D) = ("phi","q1")

#struct BoundaryConditionConstantDirichlet{T <: Real}
#    boundary_value::T
#end

#@inline function (boundary_condition::BoundaryConditionConstantDirichlet)(u, normal::AbstractVector,
#    x, t, operator_type::Trixi.Gradient,
#    equation::gravEq1D)
#return boundary_condition.boundary_value
#end

#@inline function (boundary_condition::BoundaryConditionConstantDirichlet)(u, normal::AbstractVector,
#    x, t, operator_type::Trixi.Divergence,
#    equation::gravEq1D)
#return u[2]
#end

end # module

export simpleGravity