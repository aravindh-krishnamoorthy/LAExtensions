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
function sqrtm(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("sqrt: Matrix A must be square."))
    symmetric = issymmetric(A)
    if symmetric
        e, V = eigen(A)
        negative = isreal(e) && any(e .< 0)
        if negative
            Q = Diagonal(sqrt.(complex.(e)))
            if isreal(V)
                return complex.(V*real(Q)*V', V*imag(Q)*V')
            else
                return V*Q*V'
            end
        else
            Q = Diagonal(sqrt.(e))
            return V*Q*V'
        end
    else
        U, Z = schur(A)
        # Check if lifting to complex field is necessary.
        # Lifting is only necessary if the negative diagonal element
        #   corresponds to a real eigenvalue. Negative diagonal elements
        #   corresponding to complex-conjugate pair is fine and doesn't
        #   need lifting.
        i = 1
        negative = false
        if isreal(U)
            while i < n
                if U[i,i] < 0
                    if iszero(U[i+1,i])
                        negative = true
                        break
                    else
                        i = i + 2
                        continue
                    end
                else
                    i = i + 1
                end
            end
            if i == n
                if U[i,i] < 0
                    negative = true
                end
            end
        end
        if negative
            Q = _sqrt_quasi_triu!(complex.(U))
            if isreal(Z)
                return complex.(Z*real(Q)*Z', Z*imag(Q)*Z')
            else
                return Z*Q*Z'
            end
        else
            Q = _sqrt_quasi_triu!(U)
            return Z*Q*Z'
        end
    end
end

# Square root of a 2x2 real-valued matrix with complex conjugate eigenvalues and equal diagonal values.
# Note: that T can be either <: Real or <: Complex but the values must be real, i.e., imaginary part zero.
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
function _sqrt_2x2!(A::AbstractMatrix{T}) where {T}
    @assert LinearAlgebra.checksquare(A) == 2
    @inbounds begin
        (A[1,1] == A[2,2]) || throw(ArgumentError("_sqrt_2x2!: Matrix A must have equal diagonal values."))
        p = real(A[1,2]*A[2,1])
        (p < 0) || throw(ArgumentError("_sqrt_2x2!: Matrix A must have complex conjugate eigenvalues."))
        μ = sqrt(-p)
        r = sqrt(hypot(A[1,1], μ))
        θ = atan(μ, real(A[1,1]))
        s, c = sincos(θ/2)
        α, β′ = r*c, r*s/µ
        A[1,1] = α
        A[2,2] = α
        A[1,2] = β′*A[1,2]
        A[2,1] = β′*A[2,1]
    end
    return A
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
    # Choose complex or real dot product based on T
    dot = T <: Complex ? BLAS.dotu : BLAS.dot
    # Square roots of 1x1 and 2x2 diagonal blocks
    i = 1
    sizes = ones(Int,n)
    while i < n
        if !iszero(A[i+1,i])
            _sqrt_2x2!(A[i:i+1,i:i+1])
            sizes[i] = 2
            sizes[i+1] = 0
            i += 2
        else
            A[i,i] = sqrt(A[i,i])
            sizes[i] = 1
            i += 1
        end
    end
    if i == n
        A[n,n] = sqrt(A[n,n])
        sizes[i] = 1
    end
    # Algorithm 4.3 in Reference [1]
    Δ = I(4)
    M_L₀ = zeros(T,4,4)
    M_L₁ = zeros(T,4,4)
    for k = 1:n-1
        for i = 1:n-k
            if sizes[i] == 0 || sizes[i+k] == 0 continue end
            i₁, i₂, j₁, j₂, s₁, s₂ = i, i+sizes[i]-1, i+k, i+k+sizes[i+k]-1, sizes[i], sizes[i+k]
            k₁, k₂ = i+s₁, i+k-1
            L₀ = M_L₀[1:s₁*s₂,1:s₁*s₂]
            L₁ = M_L₁[1:s₁*s₂,1:s₁*s₂]
            if s₁ == 1 && s₂ == 1
                Bᵢⱼ⁽⁰⁾ = dot(A[i₁,k₁:k₂], A[k₁:k₂,j₁])
                A[i₁,j₁] = (A[i₁,j₁] - Bᵢⱼ⁽⁰⁾)/(A[i₁,i₁] + A[j₁,j₁])
            else
                BLAS.gemm!('N', 'N', T(-1.0), A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂], T(+1.0), A[i₁:i₂,j₁:j₂])
                # Solve Uᵢ,ᵢ₊ₖ using Reference [1, (4.10)]
                kron!(L₀, Δ[1:s₂,1:s₂], A[i₁:i₂,i₁:i₂])
                L₀ .+= kron!(L₁, transpose(A[j₁:j₂,j₁:j₂]), Δ[1:s₁,1:s₁])
                ldiv!(lu!(L₀), A[i₁:i₂,j₁:j₂][:])
            end
        end
    end
    return A
end
