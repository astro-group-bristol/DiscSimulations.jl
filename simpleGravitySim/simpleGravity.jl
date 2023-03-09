module simpleGravity

using Trixi, OrdinaryDiffEq, LinearAlgebra, ForwardDiff;

struct gravEq1D <: Trixi.AbstractEquations{1 #= number of spatial dimensions =#,
                                        2 #= number of primary variables, (Φ, q1)=#}
end

@inline Trixi.flux(u, orientation, equation::gravEq1D) = SVector(0,-q1)
Trixi.varnames(::typeof(cons2cons), ::gravEq1D) = ("phi","q1")

end # module

export simpleGravity