!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2(UPLO, N, A, LDA, INFO)

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   A( LDA, * )

    EXTERNAL DTRTRS, DGEMM

    DOUBLE PRECISION   V(N)
    DOUBLE PRECISION   ONE, ZERO
    PARAMETER ( ONE = 1.0, ZERO = 0.0 )

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
        CALL DGEMM('N', 'N', I, 1, N-I, -ONE, A(1, I+1), N, A(I+1, I), N, ONE, V, I)
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
