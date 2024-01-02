using LinearAlgebra

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################

function posi!(X::AbstractMatrix{T}) where {T}
    n = size(X,1)
    v = zeros(T,n,1)
    
    ########################################
    # Cholesky decomposition
    ########################################
    LAPACK.potrf!('U', X)

    ########################################
    # MATLAB ALGORITHM
    ########################################
    # D = diag(1./diag(R)) ;
    # R(n,1:n) = conj(R(1:n,1:n)\D(1:n,n)) ;
    # for i=n-1:-1:1
    #     R(i,1:i) = conj(R(1:i,1:i)\(D(1:i,i) - R(1:i,i+1:N)*R(i+1:N,i))) ;
    # end
    # R = tril(R) ;
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
