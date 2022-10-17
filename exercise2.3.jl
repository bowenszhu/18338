using Symbolics, LinearAlgebra, FillArrays
n = 4
@variables λ[1:n]
p = prod(λ .+ 1)
s = Symbolics.scalarize(p)
expanded = expand(s)
args = arguments(Symbolics.unwrap(expanded))

@variables Λ[1:n,1:n]
Λ=Symbolics.scalarize(Λ)
u = UpperTriangular(Λ)
IΛ = u .+ Eye{Int}(n)
res = det(IΛ)
expanded = expand(res)
args = arguments(Symbolics.unwrap(expanded))
