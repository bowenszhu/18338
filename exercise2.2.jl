using Symbolics, LinearAlgebra, FillArrays
n = 4
@variables A[1:n, 1:n]
IA = A .+ Eye{Int}(n)
res = det(IA)
s = Symbolics.scalarize(res)
expanded = expand(s)
args = arguments(Symbolics.unwrap(expanded))
