!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2S(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    DOUBLE PRECISION   V(N)
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

    DO I = 1, N
        V(I) = 1/A(I,I)
    END DO
    IF (UPLO.EQ.'U') THEN
        DO J = 1, N
            A(J,J) = V(J)
            DO I = J+1, N
                A(I,J) = 0
            END DO
        END DO
        DO J = N, 1, -1
            DO I = 1, J
                DO K = N, J+1, -1
                    A(J,I) = A(J,I) - A(I,K)*A(K,J)
                END DO
            END DO
            DO K = J, 1, -1
                A(J,K) = A(J,K)*V(K)
                DO I = 1, K-1
                    A(J,I) = A(J,I) - A(I,K)*A(J,K)
                END DO
            END DO
        END DO
        DO I = 1, N
            DO J = I+1, N
                A(I,J) = A(J,I)
            END DO
        END DO
    ELSE ! UPLO.EQ.'L'
        DO J = 1, N
            A(J,J) = V(J)
            DO I = 1, J-1
                A(I,J) = 0
            END DO
        END DO
        DO J = N, 1, -1
            DO I = 1, J
                DO K = N, J+1, -1
                    A(I,J) = A(I,J) - A(K,I)*A(J,K)
                END DO
            END DO
            DO K = J, 1, -1
                A(K,J) = A(K,J)*V(K)
                DO I = 1, K-1
                    A(I,J) = A(I,J) - A(K,I)*A(K,J)
                END DO
            END DO
        END DO
        DO I = 1, N
            DO J = 1, I-1
                A(I,J) = A(J,I)
            END DO
        END DO
    END IF
    INFO = 0
    RETURN
END
