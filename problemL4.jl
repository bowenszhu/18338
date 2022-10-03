using LinearAlgebra, Distributions, Plots

n = 1000 + 1
h = 1 / (n - 1)
h² = h^2
β = 2
n_trial = 10000
largest = Vector{Float64}(undef, n_trial)
for i in 1:n_trial
    dv = randn(n)
    dv .*= 2 / sqrt(h * β)
    dv .-= range(0, 1, n)
    dv .-= 2 / h²
    ev = fill(1 / h², n - 1)
    op = SymTridiagonal(dv, ev)
    vals = eigvals(op, n:n)
    largest[i] = only(vals)
end
histogram(largest, legend = false, normalize = true)
