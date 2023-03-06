module simpleGravity

using Trixi, OrdinaryDiffEq, LinearAlgebra, ForwardDiff, Plots;

struct gravEq <: Trixi.AbstractEquations{1 #= number of spatial dimensions =#,
                                        2 #= number of primary variables, (q1, q2)=#}
end

#want to implement Φ*δ/δₜ + -q1*δ/δₓ + q2*δ/δy ₌ f(x,y)
#or Φ*δ/δₜ + -q1*δ/δₓ + q2*δ/δy - f(x,y) = 0
#where Φ = 2 + 2cos(πx)sin(2πy) and f(x,y) = 10π²cos(πx)sin(2πy).
@inline Φ(x,y) = 2 + 2cos(πx)sin(2πy)
@inline f(x,y) = 10π²cos(πx)sin(2πy)

@inline Trixi.flux(u, orientation, equation::gravEq) = ###
Trixi.varnames(_, ::gravEq) = ("scalar",)

@inline Trixi.flux_godunov(u_ll, u_rr, orientation, equation:gravEq) = flux(u_ll, orientation, equation)
@inline function Trixi.flux_ec(u_ll, u_rr, orientation, equation::gravEq)
  #return SVector(0.25 * (u_ll[1]^3 + u_ll[1]^2 * u_rr[1] + u_ll[1] * u_rr[1]^2 + u_rr[1]^3))
end

end # module

export simpleGravity