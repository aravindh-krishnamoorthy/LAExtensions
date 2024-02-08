!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2S(UPLO, N, A, LDA, INFO)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            I, J, K
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

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
        DO J = N, 1, -1
            DO CONCURRENT (I = 1:J)
                A(I,J) = A(I,J) - DOT_PRODUCT(A(J+1:N,I), A(J,J+1:N))
            END DO
            DO I = J-1, 1, -1
                A(I,J) = A(I,J) - DOT_PRODUCT(A(I+1:J,I), A(I+1:J,J))
            END DO
        END DO
        DO CONCURRENT (I = 1:N)
            A(I,1:I-1) = A(1:I-1,I)
        END DO
    END IF
    INFO = 0
    RETURN
END
