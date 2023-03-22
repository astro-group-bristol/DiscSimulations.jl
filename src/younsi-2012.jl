module PolishDoughnut

using DifferentialEquations
using StaticArrays
using ForwardDiff
using Roots

# define a derivative functor
D(f) = x -> ForwardDiff.derivative(f, x)

# metric components for the Kerr spacetime
Σ(r, θ, a) = r^2 + a^2 * cos(θ)^2
Δ(r, M, a) = r^2 - 2 * M * r + a^2
function metric_components(r, θ, M, a)
    R = 2M
    Σ₀ = Σ(r, a, θ)
    sinθ2 = sin(θ)^2

    tt = -(1 - (R * r) / Σ₀)
    rr = Σ₀ / Δ(r, R, a)
    θθ = Σ₀
    ϕϕ = sinθ2 * (r * r + a * a + (sinθ2 * R * r * a * a) / Σ₀)

    tϕ = (-R * r * a * sinθ2) / Σ₀
    @SVector [tt, rr, θθ, ϕϕ, tϕ]
end

# younsi et al. 2012 equations
# eq. (32)
Ω(z, M, a, n, rₖ) = (√M / (z^(3 / 2) + a * √M)) * (rₖ / z)^n
Ω(r, θ, M, a, n, rₖ) = Ω((r * sin(θ)), M, a, n, rₖ)

# eq. (30)
Ψ₁(r, θ, M, a, Ω, Σ) = M * ((Σ - 2r^2) / (Σ^2)) * (Ω - a * sin(θ))^2 + r * sin(θ)^2
Ψ₁(r, θ, M, a, Ω) = Ψ₁(r, θ, M, a, Ω, Σ(r, θ, a))

# eq. (31)
Ψ₂(r, θ, M, a, Ω, Σ, Δ) = sin(2θ) * ((M * r / (Σ^2)) * (a * Ω - (r^2 + a^2))^2 + Δ / 2)
Ψ₂(r, θ, M, a, Ω) = Ψ₂(r, θ, M, a, Ω, Σ(r, θ, a), Δ(r, M, a))

# eq. (28) and (29)
function drdθ(u, p, λ)
    # unpack parameters
    r, θ = u

    Ω₀ = inv(Ω(r, θ, p.M, p.a, p.n, p.rₖ))
    ψ₁ = Ψ₁(r, θ, p.M, p.a, Ω₀)
    ψ₂ = Ψ₂(r, θ, p.M, p.a, Ω₀)
    Δ₀ = Δ(r, p.M, p.a)
    Σ₀ = Σ(r, θ, p.a)

    denom = √(Δ₀ * ψ₁^2 + ψ₂^2) * √inv(Δ₀ / Σ₀)
    dr = ψ₂ / denom
    dθ = -ψ₁ / denom

    # use a static vector to avoid heap allocation
    @SVector [dr, dθ]
end

function energy(r, θ, M, a, n, rₖ)
    Ω₀ = Ω(r, θ, M, a, n, rₖ)
    # analytic solution for energy of a circular orbit
    g = metric_components(r, θ, M, a)
    -(g[1] + g[5] * Ω₀) / √(-g[1] - 2g[5] * Ω₀ - g[4] * Ω₀^2)
end

function energy_minima(θ, M, a, n, rₖ, init_r)
    # capture all but radius in a closure
    closure = r -> energy(r, θ, M, a, n, rₖ)
    # we want to find r for which dE/dr == 0
    dEdr = D(closure)
    # we give the function to find the zero for and the gradient of the function
    # to the root finder, so we can use something like the Newton solver, which
    # is pretty fast and accurate
    R = find_zero((dEdr, D(dEdr)), init_r, Roots.Newton())
    R
end

const terminate_if_negative =
    DiscreteCallback((u, t, integrator) -> u[1] * cos(u[2]) < 0, terminate!)

"""
    isobar(M, a, n, rₖ; θ = π/2, init_r = 5.0, max_λ = 40.0)

Calculate the isobar surface of a pressure supported accretion disc as in Younsi et al. (2012).
- `M`: black hole mass,
- `a`: (dimensionless) black hole spin,
- `n`: angular velocity profile index of the torus,
- `rₖ`: central Keplerian radius of the torus.

The inclination of the disc from the spin axis may be set by altering `θ`. The minima of ``E(r)`` as a function of ``r``
is numerically solved using a root finder, starting at `init_r`, and `max_λ` is the maximum integration time for the
differential equation solver.
"""
function isobar(M, a, n, rₖ; θ = π / 2, init_r = 5.0, max_λ = 40.0)
    r = energy_minima(θ, M, a, n, rₖ, init_r)
    u0 = @SVector [r, θ]

    λ_span = (0.0, max_λ)
    # named tuple syntax
    params = (; M = M, a = a, n = n, rₖ = rₖ)
    prob = ODEProblem{false}(drdθ, u0, λ_span, params)

    sol = solve(prob, Tsit5(), callback = terminate_if_negative, dtmax = 5e-2)

    # quick and dirty unpack
    r = first.(sol.u)
    θ = last.(sol.u)

    # then translate to cartesian coordinates
    z = @. r * cos(θ)

    # and only return those with z > 0
    I = @. z > 0
    r[I], z[I]
end

end # module PolishDoughnut

export PolishDoughnut
