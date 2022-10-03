using LinearAlgebra, Distributions, Plots

function T(n, β)
    dv = rand(Normal(0, √(2β)), n)
    ev = [rand(Chi(i * β)) for i in (n - 1):-1:1]
    SymTridiagonal(dv, ev)
end
function B(m, n, β)
    dv = [rand(Chi(i * β)) for i in m:-1:(m - n + 1)]
    ev = [rand(Chi(i * β)) for i in (n - 1):-1:1]
    Bidiagonal(dv, ev, :L)
end

n = 1000
β = 1.0
mc = 10000
m₁ = Vector{Float64}(undef, mc)
m₂ = Vector{Float64}(undef, mc)
for i in 1:mc
    Tₙ = T(n, β)
    vals = eigvals(Tₙ)
    m₁[i] = mean(vals)
    m₂[i] = mean(x -> x^2, vals)
end
plt₁ = histogram(m₁, legend = false, title = "the first moment")
plt₂ = histogram(m₂, legend = false, title = "the second moment")

m₃ = Vector{Float64}(undef, mc)
for i in 1:mc
    Tₙ = T(n, β)
    vals = eigvals(Tₙ)
    m₃[i] = mean(x -> x^3, vals)
end
histogram(m₃, legend = false, title = "the third moment")
