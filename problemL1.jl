using LinearAlgebra, Distributions, Statistics, Plots, LaTeXStrings

function B(n::Integer, m::Integer, β)
    @assert n ≥ m
    dv = [rand(Chi(i * β)) for i in n:-1:(n - m + 1)]
    ev = [rand(Chi(i * β)) for i in (m - 1):-1:1]
    Bidiagonal(dv, ev, :L)
end

ns = 1 .<< (5:10)
βs = exp2.(1:10)
μs = Matrix{Float64}(undef, length(ns), length(βs))
σ²s = similar(μs)
n_trial = 1000
sums = Vector{Float64}(undef, n_trial)
for (i, n) in enumerate(ns)
    for (j, β) in enumerate(βs)
        println("n = $n, β = $β")
        for k in 1:n_trial
            Bₙ = B(n, n, β)
            vals = svdvals(Bₙ)
            sums[k] = sum(vals)
        end
        μ = mean(sums)
        μs[i, j] = μ
        σ²s[i, j] = varm(sums, μ)
    end
end

pltμ = plot(xlabel = L"n", ylabel = "mean", legend = :topleft)
for j in 1:length(βs)
    plot!(pltμ, ns, view(μs, :, j), label = "β = $(βs[j])")
end
plot!(pltμ)

pltσ² = plot(xlabel = L"n", ylabel = "variance", legend = :topleft)
for j in 1:length(βs)
    plot!(pltσ², ns, view(σ²s, :, j), label = "β = $(βs[j])")
end
plot!(pltσ²)

μs = similar(ns, Float64)
for (i, n) in enumerate(ns)
    println("n = $n")
    for k in 1:n_trial
        A = randn(n, n) ./ n
        vals = svdvals(A)
        sums[k] = sum(vals)
    end
    μ = mean(sums)
    μs[i] = μ
end

pltLμ = plot(ns, μs, xlabel = L"n", ylabel = "mean", legend = false)

function mat(n::Integer, z::Complex)
    temp = randn(n, n) + im * randn(n, n)
    temp ./= sqrt(2n)
    temp[diagind(temp)] .-= z
    temp
end
n = 1000
xs = @view collect(range(0, 2, 8 + 1))[1:(end - 1)]
zs = cispi.(xs) .* 0.5
hists = Plots.Plot[]
for (i, z) in enumerate(zs)
    A = mat(n, z)
    vals = svdvals(A)
    @. vals = vals^2
    plt = histogram(vals, legend = false, title = "z = $z")
    push!(hists, plt)
end
plot(hists..., layout = (4, 2), titlefont = 6)

zs = cispi.(xs) .* 2
hists = Plots.Plot[]
for (i, z) in enumerate(zs)
    A = mat(n, z)
    vals = svdvals(A)
    @. vals = vals^2
    plt = histogram(vals, legend = false, title = "z = $z")
    push!(hists, plt)
end
plot(hists..., layout = (4, 2), titlefont = 6)
