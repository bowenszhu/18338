using LinearAlgebra, Combinatorics, Distributions, Statistics

function randprojDPP(Y)
    n = size(Y, 2)
    F = Vector{Int}(undef, n)
    for k in 1:n
        p = mean(abs2, Y, dims = 2)
        F[k] = rand(Categorical(p))
        Y = (Y * qr(Y[F[k], :]).Q)[:, 2:end]
    end
    return sort(F)
end

function randDPP(Y, Λ)
    mask = @. rand(Bernoulli(Λ / (Λ + 1)))
    return randprojDPP(Y[:, mask])
end
