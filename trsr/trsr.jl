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
function sqrtm(A::AbstractMatrix{T}, f::Function = trsr!) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("sqrt: Matrix A must be square."))
    symmetric = ishermitian(A)
    if symmetric
        e, V = eigen(A)
        negative = any(e .< 0)
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
    elseif isreal(A)
        S = schur(real(A))
        negative = false
        for i = 1:n
            if isreal(S.values[i]) && real(S.values[i]) < 0
                negative = true
                break
            end
        end
        if negative
            S = Schur{Complex}(S)
        end
        return S.Z*f(S.T)*S.Z'
    else # complex A
        S = schur(A)
        return S.Z*f(S.T)*S.Z'
    end
end

################################################################################
#
# Square root of a quasi upper triangular matrix (output of Schur decomposition)
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# VERSION: Pure Julia version for both real- and complex-valued inputs
# NOTE: It is assumed that the diagonal elements of A have a square root in type T
#
################################################################################
@views @inbounds function trsr!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("trsr!: Matrix A must be square."))
    # Choose complex or real dot product based on T
    dot = T <: Complex ? BLAS.dotu : BLAS.dot
    # Square roots of 1x1 and 2x2 diagonal blocks
    i = 1
    sizes = ones(Int,n)
    while i < n
        if !iszero(A[i+1,i])
            μ = sqrt(-real(A[i,i+1]*A[i+1,i]))
            r = sqrt(hypot(A[i,i], μ))
            θ = atan(μ, real(A[i,i]))
            s, c = sincos(θ/2)
            α, β′ = r*c, r*s/µ
            A[i,i] = α
            A[i+1,i+1] = α
            A[i,i+1] = β′*A[i,i+1]
            A[i+1,i] = β′*A[i+1,i]
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
                # Compute Bᵢⱼ⁽⁰⁾ and update A[i₁:i₂,j₁:j₂]
                mul!(A[i₁:i₂,j₁:j₂], A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂], T(-1.0), T(+1.0))
                # Solve Uᵢ,ᵢ₊ₖ using Reference [1, (4.10)]
                kron!(L₀, Δ[1:s₂,1:s₂], A[i₁:i₂,i₁:i₂])
                L₀ .+= kron!(L₁, transpose(A[j₁:j₂,j₁:j₂]), Δ[1:s₁,1:s₁])
                ldiv!(lu!(L₀), A[i₁:i₂,j₁:j₂][:])
            end
        end
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
# VERSION: Eventual FORTRAN version for both real- and complex-valued inputs
# NOTE: It is assumed that the diagonal elements of A have a square root in type T
#
################################################################################
@views @inbounds function trsrf!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("trsrf!: Matrix A must be square."))
    # Choose complex or real dot product based on T
    dot = T <: Complex ? BLAS.dotu : BLAS.dot
    # Square roots of 1x1 and 2x2 diagonal blocks
    i = 1
    sizes = ones(Int,n)
    while i < n
        if !iszero(A[i+1,i])
            μ = sqrt(-real(A[i,i+1]*A[i+1,i]))
            r = sqrt(hypot(A[i,i], μ))
            θ = atan(μ, real(A[i,i]))
            s, c = sincos(θ/2)
            α, β′ = r*c, r*s/µ
            A[i,i] = α
            A[i+1,i+1] = α
            A[i,i+1] = β′*A[i,i+1]
            A[i+1,i] = β′*A[i+1,i]
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
                # Compute Bᵢⱼ⁽⁰⁾ and update A[i₁:i₂,j₁:j₂]
                BLAS.gemm!('N', 'N', T(-1.0), A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂], T(+1.0), A[i₁:i₂,j₁:j₂])
                # Solve Uᵢ,ᵢ₊ₖ
                _, scale = LAPACK.trsyl!('N', 'N', A[i₁:i₂,i₁:i₂], A[j₁:j₂,j₁:j₂], A[i₁:i₂,j₁:j₂])
                rmul!(A[i₁:i₂,j₁:j₂], inv(scale))
            end
        end
    end
    return A
end

