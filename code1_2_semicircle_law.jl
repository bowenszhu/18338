using LinearAlgebra, Plots
n = 1000
A = randn(n, n)
Hₙ = (A + A') / 2
v = eigvals(Hₙ)
v ./= sqrt(n / 2)

dx = 0.1
x = -2:dx:2
histogram(v, bins = x, normalize = true, aspect_ratio = 2π, legend = false)
semicircle = @. sqrt(4 - x^2) / (2π)
plot!(x, semicircle, lw = 3)
savefig("semicircle.png")
