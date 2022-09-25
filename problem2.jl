using LinearAlgebra: SymTridiagonal, mul!, eigen, eigvals, Bidiagonal
using SparseArrays: sparsevec
using SpecialFunctions: gamma

function orthopoly_evaluate(T::SymTridiagonal, x::AbstractVector, isall::Bool)
    if isall
        n = size(T, 1)
        ϕ = Matrix{Float64}(undef, n, length(x))
        λ, Q = eigen(T)
        qₙ = vec(@view Q[end, :])
        temp = Vector{Float64}(undef, n)
        for (i, xᵢ) in enumerate(x)
            @. temp = λ - xᵢ
            temp .\= qₙ
            mul!(view(ϕ, :, i), Q, temp)
            ϕ[:, i] ./= ϕ[begin, i]
        end
        return ϕ
    else
        λ = eigvals(@view T[begin:(end - 1), begin:(end - 1)])
        ϕ = ones(length(x))
        for λᵢ in λ
            @. ϕ *= x - λᵢ
        end
        ϕ ./= prod(T.ev)
        return ϕ
    end
end

function compute_hermite(n::Integer, x::AbstractVector, all::Bool = true)
    ev = sqrt.(0.5:0.5:(0.5n))
    dv = zeros(n + 1)
    Tᴴ = SymTridiagonal(dv, ev)
    c = sqrt(π)
    ϕ = orthopoly_evaluate(Tᴴ, x, all)
    ϕ ./= sqrt(c)
end

function compute_laguerre(n::Integer, α::Real, x::AbstractVector, all::Bool = true)
    ev = map(i -> -sqrt(i), 1:(n - 1))
    dv = map(i -> sqrt(α + i), 1:n)
    B = Bidiagonal(dv, ev, :U)
    Tᴸ = SymTridiagonal(B' * B)
    c = gamma(α + 1)
    ϕ = orthopoly_evaluate(Tᴸ, x, all)
    ϕ ./= sqrt(c)
end

function compute_jacobi(n::Integer, α::Real, β::Real, x::AbstractVector, all::Bool = true)
    @assert α > -1 && β > -1
    dv = map(m -> sqrt((β + m) * m / ((α + β + 2m) * (1 + α + β + 2m))), 0:(n - 1))
    if α + β == -1
        dv[begin] = sqrt((α + 1) / (α + β + 2))
    end
    ev = map(m -> -sqrt((β + m) * m / ((α + β + 2m) * (1 + α + β + 2m))), 2:n)
    B = Bidiagonal(dv, ev, :U)
    Tᴶ = SymTridiagonal(B' * B)
    c = α + β == -1 : π : gamma(α + 1) * gamma(β + 1) / gamma(α + β + 1) / (α + β + 1)
    ϕ = orthopoly_evaluate(Tᴶ, x, all)
    ϕ ./= sqrt(c)
end
