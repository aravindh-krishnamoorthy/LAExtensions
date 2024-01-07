################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using LinearAlgebra

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################

function potri2!(X::AbstractMatrix{T}) where {T}
    n = size(X,1)
    v = zeros(T,n,1)   
    ########################################
    # Inversion
    ########################################
    v[n] = 1/X[n,n]
    @views LAPACK.trtrs!('U', 'N', 'N', X[1:n,1:n], v)
    X[n,1:n] = conj(v)
    for i=n-1:-1:1
        fill!(v, 0)
        v[i] = 1/X[i,i]
        @views BLAS.gemm!('N', 'N', T(-1), X[1:i,i+1:n], X[i+1:n,i], T(+1), v[1:i])
        @views LAPACK.trtrs!('U', 'N', 'N', X[1:i,1:i], v[1:i])
        X[i,1:i] = conj(v[1:i])
    end
    for i=1:n for j=i+1:n X[i,j] = X[j,i]' end end
    return X
end
