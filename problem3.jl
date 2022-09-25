using Polynomials, SpecialPolynomials

function level_density_hermite(n::Integer, xs::AbstractVector)
    xx = variable(Polynomial{Rational{Int}})
    hermite_polys = [basis(Hermite, i)(xx) for i in 0:(n - 1)]
    ϕ = Matrix{Float64}(undef, n, length(xs))
    for (i, x) in enumerate(xs)
        c = exp(-x^2 / 2)
        for (j, poly) in enumerate(hermite_polys)
            ϕ[j, i] = c * poly(x)
        end
    end
    for j in 0:(n - 1)
        ϕ[j + 1, :] .*= ((1 << j) * factorial(j) * √π)^-0.5
    end
    map(i -> sum(abs2, @view ϕ[:, i]), eachindex(xs)) # ρᴴ
end

xs = collect(-3:0.1:3)
ldᴴ = level_density_hermite(3, xs)
using Plots
plot(xs, ldᴴ)
