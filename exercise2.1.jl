using LinearAlgebra

binary2gray(n) = n ⊻ (n >> 1)

function detAB_Gray(A::AbstractMatrix, B::AbstractMatrix)
    C = copy(A)
    C_copy = copy(C) # reduce allocation due to LU
    det_sum = det(lu!(C_copy; check = false))
    gray_last = 0
    for i in 1:((1 << size(A, 2)) - 1)
        gray = binary2gray(i)
        diff = gray ⊻ gray_last
        s = trailing_zeros(diff) + 1
        if gray & diff == 0
            C[:, s] .= @view A[:, s]
        else
            C[:, s] .= @view B[:, s]
        end
        C_copy .= C
        det_sum += det(lu!(C_copy; check = false))
        gray_last = gray
    end
    det_sum
end

n = 10
A = randn(n, n)
B = randn(n, n)
@show det(A + B), detAB_Gray(A, B)
