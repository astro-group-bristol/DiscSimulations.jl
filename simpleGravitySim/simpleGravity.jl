module simpleGravity

using Trixi, OrdinaryDiffEq, LinearAlgebra, ForwardDiff;

struct gravEq1D <: Trixi.AbstractEquations{1 #= number of spatial dimensions =#,
                                        2 #= number of primary variables, (Φ, q1)=#}
end

@inline Trixi.flux(u, orientation, equation::gravEq1D) = SVector(-u[1], 0)
Trixi.varnames(::typeof(cons2cons), ::gravEq1D) = ("Φ","q1")

end # module

export simpleGravity