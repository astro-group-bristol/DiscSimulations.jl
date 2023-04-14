function central_difference(Δu, N)
    lower = -ones(Float64, N - 1)
    mid = zeros(Float64, N)
    upper = ones(Float64, N - 1)
    LinearAlgebra.Tridiagonal(lower, mid, upper) / 2Δu
end

function STANDARD_BURGER_INIT(x)
    @. 1 - cos(x)
end

const STANDARD_BURGER_TSPAN = (0.0, 1.0)