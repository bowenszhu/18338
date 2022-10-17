using Symbolics, LinearAlgebra
n = 4
@variables A[1:n, 1:n]
IA = collect(A) + I
res = det(IA)
expanded = expand(res)
args = arguments(Symbolics.unwrap(expanded))
