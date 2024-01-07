SUBROUTINE DPOTRI2(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    EXTERNAL DTRTRS, DGEMM

    DOUBLE PRECISION   V(N)
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

 !   n = size(X,1)
 !   v = zeros(T,n,1)   
 !   ########################################
 !   # Inversion
 !   ########################################
 !   v[n] = 1/X[n,n]
 !   @views LAPACK.trtrs!('U', 'N', 'N', X[1:n,1:n], v)
 !   X[n,1:n] = conj(v)
 !   for i=n-1:-1:1
 !       fill!(v, 0)
 !       v[i] = 1/X[i,i]
 !       @views BLAS.gemm!('N', 'N', T(-1), X[1:i,i+1:n], X[i+1:n,i], T(+1), v[1:i])
 !       @views LAPACK.trtrs!('U', 'N', 'N', X[1:i,1:i], v[1:i])
 !       X[i,1:i] = conj(v[1:i])
 !   end
 !   for i=1:n for j=i+1:n X[i,j] = X[j,i]' end end
 !   return X    

    DO I = 1, N
        V(I) = 0
    END DO
    V(N) = 1/A(N,N)
    CALL DTRTRS('U', 'N', 'N', N, 1, A, N, V, N, INFO)
    IF (INFO.NE.0) THEN
        RETURN
    END IF
    DO I = 1, N
        A(N,I) = V(I)
    END DO

    DO I = N-1, 1, -1
        DO J = 1, N
            V(J) = 0
        END DO
        V(I) = 1/A(I,I)
        CALL DGEMM('N', 'N', I, 1, N-I, -ONE, A(1, I+1), N, A(I+1, I), N, ONE, V, N-I)
        CALL DTRTRS('U', 'N', 'N', I, 1, A, N, V, I, INFO)
        DO J = 1, N
            A(I,J) = V(J)
        END DO   
    END DO
    DO I = 1, N
        DO J = I+1, N
            A(I,J) = A(J,I)
        END DO
    END DO

    INFO = 0
    RETURN
END
