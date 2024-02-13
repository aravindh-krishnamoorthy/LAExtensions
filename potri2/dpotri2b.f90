!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2B(UPLO, N, A, LDA, INFO)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            NB, JB
    INTEGER            I, J, K
    PARAMETER          ( NB = 32 )

    IF (UPLO.EQ.'U') THEN
        DO CONCURRENT (I = 1:N)
            A(I,I) = 1/A(I,I)
            A(I,I+1:N) = A(I,I+1:N)*A(I,I)
            A(I+1:N,I) = 0
            A(I,I) = A(I,I)*A(I,I)
        END DO
        DO J = N, 1, -1
            DO CONCURRENT (I = 1:J)
                A(J,I) = A(J,I) - DOT_PRODUCT(A(I,J+1:N), A(J+1:N,J))
            END DO
            DO I = J-1, 1, -1
                A(J,I) = A(J,I) - DOT_PRODUCT(A(I,I+1:J), A(J,I+1:J))
            END DO
        END DO
        DO CONCURRENT (I = 1:N)
            A(I,I+1:N) = A(I+1:N,I)
        END DO
    ELSE ! UPLO.EQ.'L'
        DO CONCURRENT (J = 1:N)
            A(J,J) = 1/A(J,J)
            A(J+1:N,J) = A(J+1:N,J)*A(J,J)
            A(J,J+1:N) = 0
            A(J,J) = A(J,J)*A(J,J)
        END DO
        DO J = N, 1-NB, -NB
            JB = MIN(J,NB)
            CALL DPOTRI2BDDBLK(UPLO, N, A, LDA, INFO, J, J-JB+1)           
            DO CONCURRENT (K = J-JB+1:J)
                A(1:J-JB,K) = A(1:J-JB,K) - MATMUL(A(K,K+1:N), A(K+1:N,1:J-JB))
            END DO
            DO I = J-JB, 1, -1
                DO CONCURRENT (K = J-JB+1:J)
                    A(I,K) = A(I,K) - DOT_PRODUCT(A(I+1:K,I), A(I+1:K,K))
                END DO
            END DO
        END DO
        DO CONCURRENT (I = 1:N)
            A(I,1:I-1) = A(1:I-1,I)
        END DO
    END IF
    INFO = 0
    RETURN
END

SUBROUTINE DPOTRI2BDDBLK(UPLO, N, A, LDA, INFO, BEG, END)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N, BEG, END
    DOUBLE PRECISION   A( LDA, * )
    INTEGER            I, K

    IF (UPLO.EQ.'U') THEN
    ELSE
        DO K = END, BEG, -1
            DO CONCURRENT (I = BEG:K)
                A(I,K) = A(I,K) - DOT_PRODUCT(A(K+1:N,I), A(K,K+1:N))
            END DO
            DO I = K-1, BEG, -1
                A(I,K) = A(I,K) - DOT_PRODUCT(A(I+1:K,I), A(I+1:K,K))
            END DO
        END DO
    END IF
    INFO = 0
    RETURN
END SUBROUTINE
