module simpleGravity

using Trixi, OrdinaryDiffEq, LinearAlgebra, ForwardDiff;

struct GravityEquations1D{T} <: Trixi.AbstractEquations{1,4}
    fluid::Trixi.CompressibleEulerEquations1D{T}
end
# add default constructor
GravityEquations1D(gamma) = GravityEquations1D(Trixi.CompressibleEulerEquations1D(gamma))

function Trixi.flux(u, orientation, equation::GravityEquations1D)
    # u is now (ρ, ρ_v1, ρ_e, a1)
    f = Trixi.flux(u[:3], orientation, equation.fluid)
    # add our acceleration component scaled by the relaxation time
    v1 = f[2] / f[1] + u[4]
    return SVector(f[1], v1, f[3], 0.0)
end

Trixi.varnames(::typeof(cons2cons), ::GravityEquations1D) = ("rho","rho_v1","rho_e","a1")

end # module

export simpleGravity