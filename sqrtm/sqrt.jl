################################################################################
#
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
#
################################################################################

import Base: sqrt
using LinearAlgebra

################################################################################
#
# Square root of a matrix
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# NOTE: The matrix is lifted to the complex field if the diagonal values of the
#         triangular Schur matrix does not have a square root in type T
#
################################################################################
function sqrt(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("sqrt: Matrix A must be square."))
    symmetric = issymmetric(A)
    S = schur(A)
    d = diag(S.T)
    if all(isreal(d)) && any(d .< 0)
        if issymmetric
            d = complex.(d)
        else
            S.T = complex.(S.T)
        end
    end
    if symmetric
        return S.Z * Diagonal(sqrt(d)) * S.Z'
    else
        return S.Z * _sqrt_quasi_triu!(S.T) * S.Z'
    end
end

################################################################################
#
# Square root of a quasi upper triangular matrix (output of Schur decomposition)
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# NOTE: It is assumed that the diagonal elements of A have a square root in type T
#
################################################################################
@views function _sqrt_quasi_triu!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("_sqrt_quasi_triu!: Matrix A must be square."))
    # Square roots of 1x1 and 2x2 diagonal blocks
    i = 1
    sizes = ones(Int,n)
    while i < n
        if !iszero(A[i+1,i])
            LinearAlgebra._sqrt_real_2x2!(A[i:i+1,i:i+1], A[i:i+1,i:i+1])
            sizes[i] = 2
            sizes[i+1] = 0
            i += 2
        else
            A[i,i] = sqrt(A[i,i])
            i += 1
        end
    end
    if i == n
        A[n,n] = sqrt(A[n,n])
    end
    # Algorithm 4.3 in Reference [1]
    Δ = I(4)
    M_L₀ = zeros(T,4,4)
    M_L₁ = zeros(T,4,4)
    M_Bᵢⱼ⁽⁰⁾ = zeros(T,2,2)
    for k = 1:n-1
        for i = 1:n-k
            if sizes[i] == 0 || sizes[i+k] == 0 continue end
            k₁, k₂ = i+1+(sizes[i+1]==0), i+k-1
            i₁, i₂, j₁, j₂, s₁, s₂ = i, i+sizes[i]-1, i+k, i+k+sizes[i+k]-1, sizes[i], sizes[i+k]
            L₀ = M_L₀[1:s₁*s₂,1:s₁*s₂]
            L₁ = M_L₁[1:s₁*s₂,1:s₁*s₂]
            Bᵢⱼ⁽⁰⁾ = M_Bᵢⱼ⁽⁰⁾[1:s₁, 1:s₂]
            # Compute Bᵢⱼ⁽⁰⁾
            mul!(Bᵢⱼ⁽⁰⁾, A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂])
            # Solve Uᵢ,ᵢ₊ₖ using Reference [1, (4.10)]
            kron!(L₀, Δ[1:s₂,1:s₂], A[i₁:i₂,i₁:i₂])
            L₀ .+= kron!(L₁, transpose(A[j₁:j₂,j₁:j₂]), Δ[1:s₁,1:s₁])
            A[i₁:i₂,j₁:j₂] .-= Bᵢⱼ⁽⁰⁾
            ldiv!(lu!(L₀), A[i₁:i₂,j₁:j₂][:])
        end
    end
    return A
end
