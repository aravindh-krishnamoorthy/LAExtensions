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

    INTEGER            I, J, K, L
    INTEGER            NN, JB, NB
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )
    PARAMETER ( NB = 32 )

    IF (UPLO.EQ.'U') THEN
        DO I = 1,N
            A(I,I) = 1/A(I,I)
            DO CONCURRENT (J = I+1:N)
                A(I,J) = A(I,J)*A(I,I)
                A(J,I) = 0
            END DO
            A(I,I) = A(I,I)*A(I,I)
        END DO
        NN = ((N-1)/NB)*NB+1
        DO J = NN, 1, -NB
            JB = MIN(NB, N-J+1)
            DO L = J+JB-1, J, -1
                DO K = N, L+1, -1
                    DO CONCURRENT (I = J:L)
                        A(L,I) = A(L,I) - A(I,K)*A(K,L)
                    END DO
                END DO
                !A(L,J:L) = A(L,J:L) - MATMUL(A(L+1:N,L), TRANSPOSE(A(J:L,L+1:N)))
                DO K = L, J, -1
                    DO CONCURRENT (I = 1:K-1)
                        A(L,I) = A(L,I) - A(I,K)*A(L,K)
                    END DO
                    !A(L,1:K-1) = A(L,1:K-1) - A(1:K-1,K)*A(L,K)
                END DO
            END DO
            DO L = J, J+JB-1
                !A(L,1:J-1) = A(L,1:J-1) - MATMUL(A(L+1:N,L), TRANSPOSE(A(1:J-1,L+1:N)))
                DO I = 1, J-1
                    DO K = L+1, N
                        A(L,I) = A(L,I) - A(I,K)*A(K,L)
                    END DO
                END DO
            END DO
            DO L = J, J+JB-1
                DO I = J-2, 1, -1
                    DO K = I+1, J-1
                        A(L,I) = A(L,I) - A(I,K)*A(L,K)
                    END DO
                END DO
            END DO
        END DO
        DO CONCURRENT (I = 1:N)
            DO CONCURRENT (J = I+1:N)
                A(I,J) = A(J,I)
            END DO
        END DO
    ELSE ! UPLO.EQ.'L'
        DO J = 1, N
            A(J,J) = 1/A(J,J)
            DO CONCURRENT (I = J+1:N)
                A(I,J) = A(I,J)*A(J,J)
                A(J,I) = 0
            END DO
            A(J,J) = A(J,J)*A(J,J)
        END DO
        NN = ((N-1)/NB)*NB+1
        DO J = NN, 1, -NB
            JB = MIN(NB, N-J+1)
            DO L = J+JB-1, J, -1
                DO K = N, L+1, -1
                    DO CONCURRENT (I = J:L)
                        A(I,L) = A(I,L) - A(K,I)*A(L,K)
                    END DO
                END DO
                DO K = L, J, -1
                    DO CONCURRENT (I = 1:K-1)
                        A(I,L) = A(I,L) - A(K,I)*A(K,L)
                    END DO
                END DO
            END DO
            DO L = J+JB-1, J, -1
                DO K = N, L+1, -1
                    DO CONCURRENT (I = 1:J-1)
                        A(I,L) = A(I,L) - A(K,I)*A(L,K)
                    END DO
                END DO
            END DO
            DO K = J-1, 1, -1
                DO I = 1,K-1
                    DO L = J+JB-1, J, -1
                        A(I,L) = A(I,L) - A(K,I)*A(K,L)
                    END DO
                END DO
            END DO
        END DO
        DO CONCURRENT (I = 1:N)
            DO CONCURRENT (J = 1:I-1)
                A(I,J) = A(J,I)
            END DO
        END DO
    END IF
    INFO = 0
    RETURN
END
